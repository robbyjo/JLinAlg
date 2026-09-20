# Replicated single-cell RNA analysis

The `single-cell` command connects quantified, annotated RNA counts to QC,
pseudobulk, differential state, relative population abundance, gene-set scores,
and exploratory PCA/clustering. Biological samples, not cells, determine the
replication used for condition inference. This is a source-build addition;
build `./gradlew check executableJar` (Windows: `.\gradlew.bat check executableJar`).

These operations implement the replicated, annotated-input RNA branch of the
roadmap. The complete protein, ATAC, trajectory, reference-annotation and
perturbation suites remain open in [TODO.md](../../TODO.md). See the
[validation report](../cell-spatial-validation.md) for exactly what was tested.

## Four explicit input tables

The runnable synthetic fixture in `examples/cell-spatial` has eight samples,
320 cells, 31 genes, two annotated populations and two tissue sections per sample.
Its values are artificial and demonstrate file contracts and effect recovery.

| Input | Required columns | Contract |
|---|---|---|
| `--counts` | `obs_id,feature_id,count` | TSV in exactly this order; nonnegative integer raw RNA counts, unique cell/gene pairs. Omitted entries are measured zeros. |
| `--cells` | `obs_id,sample_id,cell_type` | Unique observation IDs, exactly matched samples; explicit `unknown` labels are allowed. Optional `exclude` must be `true` or `false`. |
| `--samples` | `sample_id,donor_id,condition,batch,source,assay,visit` | Unique sample IDs. Assay is `scRNA` or `spatialRNA`. Numeric covariates may be added. |
| `--features` | `feature_id,mitochondrial` | Complete measured feature universe; mitochondrial flag is `true` or `false`. |

Counts are read as sparse entries without expanding the cell-by-gene matrix.
All metadata columns are preserved, including supplied annotation provenance,
reference version, confidence, segmentation quality and collection information.
Preservation does not validate those annotations. Blank labels, duplicate IDs,
unknown count IDs, fractional counts, negative counts and nonfinite values fail.
Do not encode missing/unmeasured proteins as sparse zeros or use this RNA branch
for antibody, mass-spectrometry, cytometry or accessibility measurements.

For real studies, choose sample IDs at the biological sample/visit level. Record
donor IDs even when several specimens, sections or visits have different sample
IDs. A repeated donor cannot enter an independent-sample condition comparison.

## QC and aggregation tutorial

Every output path must be a fresh directory. Run from the repository root:

```shell
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method qc --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --min-features 5 --max-mito 0.2 --out build/cell-demo/qc
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method pseudobulk --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --min-cells 10 --out build/cell-demo/pseudobulk
```

`cell-qc.tsv` records library counts, detected features, mitochondrial fraction,
retention and every applied exclusion reason. `sample-qc.tsv` reports retained
observations and donor identity. Zero libraries always fail QC. Defaults are
`--min-counts 1 --min-features 1 --max-mito 1`; these permissive defaults are not
a universal biological QC recommendation. Each operation reapplies its own stated
QC thresholds; use identical arguments or project configuration for comparable runs.

Empty droplets, ambient RNA and doublets require appropriate upstream assessment.
Import the reviewed decision through `exclude` and preserve method/version/score
columns. The command does not estimate those quantities or silently call another
engine. Inspect donor representation and rare populations before inference.

Aggregation sums retained raw counts within `(sample_id,cell_type)`; all sections
of one sample contribute to one aggregate. `pseudobulk-counts.tsv` is feature-major,
and `pseudobulk-samples.tsv` supplies the aggregate-to-sample/donor mapping, cell
counts, library counts and eligibility. No pseudocount is added to this export.
`--group-column domain` uses a supplied `domain` column instead of cell type.
Counts are conserved over retained observations. Absent populations are absent
aggregates; they never become zero-expression artificial replicates.

## Differential state: independent or paired donors

Install R and limma separately in a chosen library, for example using
`BiocManager::install("limma")` in R. The CLI does not install packages or contact
the network. Supply the executable and library paths appropriate to your host:

```shell
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method state --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --reference control --tested case --rscript Rscript --out build/cell-demo/state
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method state --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/paired-samples.tsv --features examples/cell-spatial/features.tsv --reference control --tested case --paired true --rscript Rscript --out build/cell-demo/paired-state
```

On the validated Windows host, the executable was
`"C:/Program Files/R/R-4.6.1/bin/Rscript.exe"` and the optional library argument
was `--r-library E:/Projects/JLinAlg/build/r-library`. These are host-specific paths.

Each population is modeled separately. The design includes an intercept, tested
minus reference indicator, optional numeric `--covariates age,sex_numeric`, batch
indicators, and donor fixed effects for pairs. At least two donors per group or
two complete pairs are required, with positive residual degrees of freedom.
Two donors is an identifiability minimum, not a power recommendation. Batch
effects constant within a donor pair are absorbed by donor effects. A batch that
changes within pairs is modeled explicitly. Redundant or confounded designs fail
that population with a recorded status. Source and assay must be identical within
each modeled population. Exactly two declared conditions are supported.

For paired designs, each retained donor must have exactly one sample per condition.
Incomplete pairs are reported as unsupported; the tool does not silently discard
unmatched samples. Arbitrary longitudinal data, multiple same-condition visits,
multi-level random effects and cell-level hierarchical inference remain unsupported.
Between-donor covariates that are redundant with paired donor effects must not be
included as separate columns.

The RNA engine is **limma voom**, followed by `lmFit` and `eBayes`. It uses total
assayed-library normalization, retaining the full library size after feature
filtering. No TMM or composition-factor normalization is applied. Global expression
or composition shifts can therefore affect the estimand; assess that assumption
before applying this path. Genes require count >=10 in at least two eligible
samples by default (`--min-gene-count`, `--min-gene-samples`). At least two retained
genes are needed for voom; small panels give limited mean/variance information.

`results.tsv` contains every feature/population hypothesis, log2-scale effect,
SE, 95% interval, moderated degrees of freedom, statistic, p-value, BH value and
status. BH spans populations and genes together. Filtered or unavailable tests
enter the family conservatively as p=1, while their public inference columns stay
`NA`. Entire populations with zero retained cells appear only in the input/QC audit.
BH assumes its usual independence/positive-dependence conditions; it is not a
guarantee under arbitrary cross-gene dependence.

Each `group-N` directory records its actual design, counts, full library sizes,
log2 CPM values, precision weights, R session, adapter source and external log.
`population-status.tsv` maps labels to tested or unsupported populations. Failed
R execution aborts publication; it never falls back to a different estimator.

## Relative population abundance and functional scores

```shell
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method abundance --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --reference control --tested case --out build/cell-demo/abundance
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method pathways --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --reference control --tested case --gene-sets examples/cell-spatial/gene-sets.tsv --out build/cell-demo/pathways
```

Abundance computes each sample's population count divided by **all retained cells**,
including `unknown` labels. Missing labels in an otherwise observed sample have
zero population counts. It fits sample-level Gaussian models of
`asin(sqrt(proportion))` with the same design rules; samples receive equal weight.
The effect is on this transformed relative-composition scale, not a log odds ratio
or absolute cell abundance. `proportions.tsv` retains numerator and denominator.
This is not Milo or a multinomial/Dirichlet model. Few cells, rare populations,
heterogeneous sampling precision and boundary-heavy proportions need sensitivity
analysis; a Gaussian error model may be unsuitable.

Pathway inputs have unique `gene_set,feature_id` memberships in the same namespace
as the measured panel. A score is the mean of
`log2((pseudobulk_count+0.5)*1e6/(full_library+1))` over at least two measured members.
The sample-level score is tested with Gaussian regression, preserving donor
replication and within-set gene correlation through aggregation. Missing panel
members are reported. Scores represent prespecified mean expression, not inferred
regulator activity, GSEA, competitive enrichment or pathway activation probabilities.
Overlapping gene sets are dependent. Use fixed, independently specified sets.

For hit-based [enrichment](enrichment.md), select one population and contrast,
use `status == 'tested'` as eligibility, retain all eligible nonsignificant genes
in the measured background, and select hits using the reported BH column. For
[replication/meta-analysis](meta-analysis.md), harmonize population definitions,
contrast direction, units and measured features, then export `feature_id`, `beta`
(the state effect) and `se` separately for each independent cohort. Cohort overlap,
different tissues or unmatched assay scales require additional models. Do not pool
cell numbers or sections as if they were independent studies.

## Exploratory representation

```shell
java -jar build/cli/jlinalg-0.3.6.jar single-cell --method representation --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --rscript Rscript --variable-features 20 --components 6 --clusters 3 --seed 1 --out build/cell-demo/representation
```

This branch uses `log1p(10000*count/library)`, selects features by sample variance
across cells, centers without scaling, computes PCA, and fits seeded k-means with
20 starts. It exports normalized values, selected features, loadings, PCA variance,
embedding and cluster assignments. It does not estimate batch integration, UMAP,
nearest-neighbor clustering, reference labels or marker p-values. Check batch and
condition representation before interpreting clusters. Never use these embedding
coordinates as raw counts in the differential-state command. Analyses selected
using the same clusters/features are exploratory.

## Provenance and resource limits

Project-over-local configuration is shared with every other command. Each result
directory is staged and published only on success, with `manifest.yaml` containing
input/output SHA-256 hashes, effective options and declared assumptions;
`OUT.config.yaml` records resolved configuration. Original metadata and feature
universes are copied into the result. R versions are recorded per adapter run.

Limits: 100,000 observations; 50,000 genes; 10,000 samples; two million sparse
entries; 32 MB per metadata table; ten million dense aggregate entries. Exploratory
PCA is limited to two million dense cell-by-feature entries. These are enforced
bounds, not demonstrated maximum-scale benchmarks. Sparse TSV interchange is
supported; native H5AD/MuData/QFeatures import, resumability and distributed execution
remain future work.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Helena L. Crowell et al. (2020) — muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data](../CITATIONS.md#crowell-muscat-2020) — [PMID: 33257685](https://pubmed.ncbi.nlm.nih.gov/33257685/)
- [Jordan W. Squair et al. (2021) — Confronting false discoveries in single-cell differential expression](../CITATIONS.md#squair-pseudoreplication-2021) — [PMID: 34584091](https://pubmed.ncbi.nlm.nih.gov/34584091/)
- [Gordon K. Smyth (2004) — Linear models and empirical Bayes methods for assessing differential expression in microarray experiments](../CITATIONS.md#smyth-limma-2004) — [PMID: 16646809](https://pubmed.ncbi.nlm.nih.gov/16646809/)
- [Charity W. Law, Yunshun Chen, Wei Shi, and Gordon K. Smyth (2014) — voom: precision weights unlock linear model analysis tools for RNA-seq read counts](../CITATIONS.md#law-voom-2014) — [PMID: 24485249](https://pubmed.ncbi.nlm.nih.gov/24485249/) · [PMCID: PMC4053721](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053721/)
- [Yoav Benjamini and Yosef Hochberg (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](../CITATIONS.md#benjamini-hochberg-1995)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
