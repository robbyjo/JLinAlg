# Networks, modules and candidate regulators

Six CLI methods connect completed omics results to reference networks, molecular
modules, conditional associations, regulatory predictions, signaling hypotheses
and differential connectivity. They do not label inferred hubs as proven causal
drivers. These are source-build additions; build `./gradlew check executableJar`.
Every `--out` is a **new directory**. No LLM is involved.

Every method exports `network.graphml` for Cytoscape/GraphML viewers as well as
TSV tables. Edge attributes retain their names and meanings; matrix-based graphs
include measured features even when no selected/exported edge touches them.
GraphML attributes are strings so missing inference values remain explicit.

## Choose the input that answers the question

| Method | Inputs | Engine and interpretation |
|---|---|---|
| `reference` | Typed edge table, selected genes | Native weighted propagation and degree-matched candidate-driver enrichment |
| `wgcna` | Normalized molecular matrix; optional phenotype and replication | WGCNA signed coexpression modules, eigengenes, membership/connectivity and preservation |
| `sparse` | Normalized molecular matrix | Native Gaussian neighborhood LASSO; conditional associations |
| `regulatory` | Expression matrix and candidate transcription factors | GENIE3; directed predictive feature-importance edges, unknown activation/inhibition |
| `signaling` | Signed directed prior and regulator activities | Inverse CARNIVAL; networks compatible with measured regulator activity |
| `differential` | Two independent cohorts with the same features | Native Pearson-correlation differences and Fisher-z inference |

Results-only inputs support `reference`; appropriately derived regulator activity
summaries can support `signaling`. Association betas and p-values cannot reconstruct
a coexpression network. Matrix inputs have `sample_id` first, then numerical
feature columns, one independent biological sample per row. This orientation is
explicit and can differ from upstream JLinAlg feature-by-sample outputs: transpose
those before use, preserving IDs. All features must be finite and variable.

Normalize and handle batch/covariates before network construction, preserving
the biological contrast of interest. Raw RNA-seq counts are not accepted as a
statistically appropriate substitute for transformed expression. Do not select
only differential hits before WGCNA; overlay discovery results on modules built
from an appropriately filtered measured universe. Repeated visits/cells from one
donor do not become independent samples by giving them new IDs. Single-cell
donor inference needs an appropriate aggregation/design outside these commands.

The current matrix implementation accepts at most 2,000 features and defaults to
10,000 samples, with an additional 20-million-cell bound. `--max-features` can
lower that limit; `--max-samples` can raise the sample limit within the cell bound.
These are bounded follow-up workflows, not a validated 30,000-gene dense engine.
Filter by an outcome-blind rule and document the feature universe.

## Reference networks and candidate key drivers

```shell
java -jar build/cli/jlinalg-0.3.6.jar network --method reference --edges examples/followup/edges.tsv --hits examples/followup/hits.tsv --permutations 999 --seed 19 --out build/followup/reference
```

Edges require `source,target,weight,type`, with positive weights, distinct nodes,
and no duplicate pair. Delimit fields with tabs. Set `--directed true` for directed
graphs; otherwise edges are undirected. Choose/aggregate an evidence layer
explicitly: duplicate evidence records are not automatically summed. Extra input
columns survive in `edges.tsv`. STRING functional associations are not all direct
physical interactions; select and label the desired resource/network layer.

Hits require a unique `gene` column using the same namespace. `unmapped.tsv`
records selected genes absent from the network; they are not treated as negatives.
The network universe consists of nodes appearing in edges; isolated measured genes
are not automatically added. At least one selected gene must map.

`nodes.tsv` reports degree, seed status, propagation rank and candidate-driver
tests. Propagation is a weighted random walk with restart (default 0.5), starting
from equal mass on mapped hits. Dangling-node mass returns to that seed distribution.
It uses edge magnitudes and does **not** propagate activation/inhibition signs.

Candidate-driver enrichment counts selected genes in each node's outgoing
neighborhood (undirected neighborhood by default), excluding the candidate itself.
`--radius` is 1 by default and at most 3. Seed labels are permuted within degree
bins: degree 0, degree 1, degrees 2–3, 4–7, etc. Each permutation preserves the
number of seeds in every bin. The empirical one-sided p-value is `(exceedances+1)/
(permutations+1)`. BH covers **all** network nodes; overlapping neighborhoods can
induce dependence, so BH is not a blanket guarantee under arbitrary dependence.
`--permutations 0` disables tests and reports NA. Small bins can make the null
uninformative; this is a documented coarse degree-matched null, not an exact test
of biological causation. The minimum p-value is 1/(B+1).

The graph is limited to 100,000 nodes/1,000,000 edges, and neighborhood permutation
work to 100 million visits. Propagation must converge within `--iterations` at the
specified `--tolerance`; otherwise the run fails.

## WGCNA modules, hubs, phenotype association and preservation

Install R packages separately, into a chosen library if desired:

```r
install.packages("BiocManager")
BiocManager::install(c("WGCNA", "GENIE3", "CARNIVAL"), ask=FALSE, update=FALSE)
```

Specify the Rscript executable through CLI or local configuration. On Windows,
use a quoted absolute path such as `C:/Program Files/R/R-4.6.1/bin/Rscript.exe`.
`--r-library DIRECTORY` prepends a separately installed library; no package
installation or source download happens during analysis.

```shell
java -jar build/cli/jlinalg-0.3.6.jar network --method wgcna --matrix examples/followup/discovery.tsv --power 6 --min-module-size 5 --pheno examples/followup/phenotype.tsv --trait trait --reference-matrix examples/followup/replication.tsv --permutations 100 --rscript Rscript --out build/followup/wgcna
```

This calls WGCNA `blockwiseModules` with signed Pearson adjacency and signed TOM,
one block and one thread. Soft-threshold power is **required**, not selected to
maximize disease association. Inspect network diagnostics and justify it for your
data. The example uses small synthetic modules; the default minimum module size
is 30, and merge cut height is 0.25. At least 15 independent samples are required;
20+ is recommended and neither threshold guarantees adequate power.

Outputs include `modules.tsv` (module, signed eigengene membership, within-module
adjacency connectivity), `eigengenes.tsv`, and `edges.tsv` (adjacency above
`--edge-threshold`, default 0.1). The edge threshold changes the exported graph,
not module discovery. Grey/unassigned genes have no within-module hub score.
Eigengene signs are arbitrary: retain them consistently within a run.

Optional phenotype input uses `sample_id` and a complete variable numeric trait.
It is aligned to matrix rows; `module-trait.tsv` contains unadjusted Pearson
associations, Student-t p-values and BH across tested eigengenes. Covariate/cluster
adjustment must be handled upstream; selecting modules/features using the same
phenotype makes subsequent p-values exploratory.

An independent `--reference-matrix` must have exactly the same feature set;
columns are reordered by ID and overlapping sample IDs are rejected. WGCNA
`modulePreservation` writes `preservation.tsv` plus the full `preservation.rds`.
Its `gold` row is WGCNA's random reference module, not a discovered biological
module. Preservation Z summaries are not association p-values. Permutation count
and random seed are explicit; the small test fixture uses 10 permutations only
to exercise the interface, not as an adequate scientific analysis.

## Sparse conditional-association networks

```shell
java -jar build/cli/jlinalg-0.3.6.jar network --method sparse --matrix examples/followup/discovery.tsv --lambda 0.15 --rule and --bootstraps 20 --seed 19 --out build/followup/sparse
```

Every feature is centered/scaled to population variance one and regressed on the
remaining features with Gaussian LASSO. Positive lambda is required and should be
chosen without optimizing a desired biological finding. `and` retains an edge
when both neighborhood fits select it; `or` requires either fit. `edges.tsv`
retains both directed regression coefficients, their mean absolute magnitude
and optional bootstrap selection frequency. Coefficients are not partial
correlations; an edge is undirected conditional-association evidence.

No p-values are fabricated after selection. Bootstrap frequency is stability,
not posterior probability. A nonconverged fit or degenerate bootstrap matrix
fails the run rather than dropping inconvenient replicates. Work is bounded by
500 million sample-feature-pair visits over the requested fits; that count does
not include coordinate-descent iterations and is not a runtime guarantee.
The public Java API is `org.jlinalg.network.NetworkAnalysis.neighborhoodLasso`.

## Regulatory networks with GENIE3

```shell
java -jar build/cli/jlinalg-0.3.6.jar network --method regulatory --matrix examples/followup/discovery.tsv --regulators examples/followup/regulators.tsv --trees 1000 --seed 19 --rscript Rscript --out build/followup/regulatory
```

The regulator file has a unique `gene` column; every candidate regulator must be
measured. The adapter transposes to GENIE3's feature-by-sample orientation and
uses one core for seeded reproducibility. Weights are predictive importances,
not p-values or causality probabilities. Edge signs remain `unknown`; importance
does not infer activation/inhibition. `weights.rds` preserves the full result.
This adapter is GENIE3, not a claim of SCENIC motif pruning, SCENIC+ multiome
support, or GRNBoost2 equivalence. Those are separate future adapters.

## Contextual signaling with inverse CARNIVAL

```shell
java -jar build/cli/jlinalg-0.3.6.jar network --method signaling --edges examples/followup/signaling.tsv --activities examples/followup/activities.tsv --solver lpSolve --rscript Rscript --out build/followup/signaling
```

The prior has `source,target,sign`, sign +1/-1, unique directed pairs, and no self
loops. The activity table has `gene,activity`, with discretized regulator activity
+1/-1; it must not substitute raw expression fold changes for regulatory activity
without justification. Every measured regulator must appear in the prior.

The adapter calls `runInverseCarnival`, suitable when perturbation targets are
unknown. lpSolve is useful for small examples. For larger networks, supply
`--solver cbc` or `cplex` and `--solver-path EXECUTABLE`; users install those solvers
separately. The current validation used lpSolve; CBC/CPLEX execution is not
established by that test. `--timeout` bounds the external process.

`edges.tsv` contains selected biological edges with canonical source/target/sign,
solution frequency divided by 100, and type `contextual_signaling`. Frequencies
are across returned solver solutions, not calibrated causal probabilities.
Artificial `Perturbation` edges are separated into `proposed-perturbations.tsv`.
`weightedSIF.tsv`, `nodes.tsv`, `solution.rds`, and solver files preserve the
unfiltered engine evidence. An inferred upstream regulator is a hypothesis for
perturbation follow-up, not an experimentally established driver.

## Differential connectivity

```shell
java -jar build/cli/jlinalg-0.3.6.jar network --method differential --matrix examples/followup/discovery.tsv --matrix-b examples/followup/replication.tsv --out build/followup/differential
```

This example compares independent synthetic cohorts. A biological application
would use appropriately matched conditions/cohorts. Matrices must contain the
same feature set and disjoint sample IDs; at least four independent samples per
group are required, with larger samples preferable for the approximation.

For each feature pair, z = `(atanh(r_a)-atanh(r_b))/sqrt(1/(n_a-3)+1/(n_b-3))`.
Two-sided normal p-values assume approximately bivariate-normal independent
observations. `edges.tsv` contains both correlations, their difference, z, p and
BH. Perfect/nearly perfect correlations are marked `degenerate_correlation`
with NA inference and count at p=1 in the complete pairwise correction family.
This is differential **marginal correlation**, not differential precision-matrix
inference, paired analysis or a test of causal rewiring.

## Reproducible validation and publication

`tools/generate_followup_examples.py` regenerates the seeded synthetic examples.
`tools/generate_network_reference.R` generates independent glmnet neighborhood
coefficients and base-R Fisher-z/BH fixtures in `src/test/resources/network-reference`.
`FollowupWorkflowTest` checks native methods against them, plus analytic LASSO,
propagation conservation, input alignment and failure boundaries.

All three actual R adapters were run on the synthetic examples using WGCNA 1.74,
GENIE3 1.34.0 and CARNIVAL 2.22.0 with lpSolve. Session files record exact package
versions. This validates the exercised adapters, not biological accuracy or
performance on whole-transcriptome datasets. Run `tools/validate_followup.py` with
the executable JAR, Rscript and optional R library to repeat the adapter checks.

For a paper, archive manifests, configuration sidecars, node/edge tables and the
R session information. Gene/module/variant scores remain separate evidence
dimensions. Networks can feed [enrichment](enrichment.md) or target follow-up, but
they do not replace [MR](mendelian-randomization.md), [colocalization](colocalization.md),
functional validation or target-safety assessment.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Peter Langfelder and Steve Horvath (2008) — WGCNA: an R package for weighted correlation network analysis](../CITATIONS.md#langfelder-wgcna-2008)
- [Peter Langfelder et al. (2011) — Is My Network Module Preserved and Reproducible?](../CITATIONS.md#langfelder-preservation-2011)
- [Nicolai Meinshausen and Peter Buhlmann (2006) — High-dimensional graphs and variable selection with the Lasso](../CITATIONS.md#meinshausen-networks-2006)
- [Vân Anh Huynh-Thu et al. (2010) — Inferring Regulatory Networks from Expression Data Using Tree-Based Methods](../CITATIONS.md#huynh-thu-genie3-2010)
- [Anika Liu et al. (2019) — From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL](../CITATIONS.md#liu-carnival-2019)
- [R. A. Fisher (1921) — On the probable error of a coefficient of correlation deduced from a small sample](../CITATIONS.md#fisher-correlation-1921)
- [Sergey Brin and Lawrence Page (1998) — The anatomy of a large-scale hypertextual Web search engine](../CITATIONS.md#brin-pagerank-1998)
- [Yoav Benjamini and Yosef Hochberg (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](../CITATIONS.md#benjamini-hochberg-1995)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
