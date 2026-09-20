# Physical spatial RNA analysis

The `spatial` command operates on quantified RNA measurements and supplied physical
coordinates. It builds graphs, tests within-section spatial patterns and label
relationships, compares sample-level spatial composition across conditions, and
compares distance gradients. This is distinct from regional EWAS, whose distance
is along a chromosome. See [single-cell inputs](single-cell.md) and the
[validation report](../cell-spatial-validation.md).

## Input, geometry and QC

Use the same four count/cell/sample/feature tables as `single-cell`. Add these
columns to the cell table:

| Column | Meaning |
|---|---|
| `specimen_id` | Belongs to exactly one biological sample. |
| `section_id` | Globally unique; belongs to exactly one specimen. |
| `compartment` | Supplied anatomical boundary group; graphs do not cross it. |
| `stratum` | Exchangeability group for label permutations within compartment. |
| `x,y` and optional `z` | Finite physical coordinates in one declared coordinate system. |

Specify `--units micrometer` or `millimeter`, `--coordinate-system` with a versioned
coordinate-frame description, and `--resolution cell` or `spot`. Coordinates must
already have their intended transformation applied. Preserve image references,
segmentation IDs/quality, mask provenance and transform metadata in extra columns.
QC exclusions use the same RNA metrics and optional imported `exclude` flag as
the single-cell command. Subcellular molecules must be assigned/quantified upstream.

The implementation does not register images, segment cells, estimate transforms,
read pixel masks or infer cell mixtures. Separate sections remain disconnected
even when their numeric coordinates overlap or when z values suggest adjacency.
To respect irregular holes, use a reviewed supplied adjacency graph or explicitly
partition the compartment labels. A radius alone cannot detect a hole crossed by
a short straight edge.

## Build graphs and inspect tissue maps

```shell
java -jar build/cli/jlinalg-0.3.6.jar spatial --method graph --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --units micrometer --coordinate-system synthetic-v1 --resolution cell --graph radius --radius 10.1 --out build/spatial-demo/graph
```

Radius edges connect observations at Euclidean distance <= the radius within the
same section and compartment. Two- and three-dimensional distances are supported.
Use `--graph knn --k 6 --radius 30` for an undirected **union** kNN graph with a
hard maximum gap: an edge exists if either endpoint selects the other. Ties follow
input row order. Graph construction is deterministic for the same input order.

Alternatively use `--graph adjacency --edges reviewed-edges.tsv`, omitting radius
and k. The edge table requires `source,target` observation IDs; duplicate undirected
pairs, self-edges and cross-compartment/section edges fail. Nodes excluded by QC
are removed, with their exclusions retained in `cell-qc.tsv`. Supplied edges are
unweighted; measured distances are exported for auditing.

`edges.tsv` contains each unordered pair once. `nodes.tsv` records degree and
isolates. `maps.tsv` maps each section ID to a separate SVG; circles include cell
ID/type tooltips, and isolated nodes are red. Three-dimensional inputs are shown
as explicitly labeled XY projections. Axes preserve aspect ratio; each map states
its original coordinate range and physical units. Compare radius/k choices and
compartment annotations before interpreting a graph.

## Spatially variable features

```shell
java -jar build/cli/jlinalg-0.3.6.jar spatial --method autocorrelation --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --units micrometer --coordinate-system synthetic-v1 --resolution cell --radius 10.1 --feature-list examples/cell-spatial/selected-features.tsv --permutations 999 --seed 19 --out build/spatial-demo/autocorrelation
```

For each section and requested feature, the response is
`log1p(10000*count/library)`. Moran's I and Geary's C use symmetric binary graph
weights. With n section observations, m undirected edges and centered values z:

```text
I = n * sum_edges(z_i*z_j) / (m * sum_i(z_i^2))
C = (n-1) * sum_edges((y_i-y_j)^2) / (2*m * sum_i(z_i^2))
```

Isolates remain in the section universe and denominator. Fewer than three
observations, no edges or constant values give an explicit unavailable status.
Centered values are scaled before products to avoid overflow/underflow; constant
decimal measurements are detected before rounding in the mean can mimic variation.
`--feature-list` requires a `feature_id` column; otherwise all measured genes are
requested. Prefer a prespecified feature family when permutations would be large.

Random-label permutations operate within section/compartment/stratum, holding
geometry fixed. Tests use the fixed two-sided statistics `abs(I+1/(n-1))` and
`abs(C-1)`. Those centers are unrestricted random-label expectations; restricted
null distributions can have different means, but permutation validity uses the
same statistic for observed and permuted data. Exchangeability must hold within
each supplied stratum. A specimen with singleton strata has no permutation power.

The reported Monte Carlo p-value is `(exceedances+1)/(B+1)`, with minimum resolution
`1/(B+1)`. BH covers both statistics, all requested genes and all sections together;
unavailable rows contribute p=1 internally. This tests within-section pattern,
not population-level differential expression or a condition effect. It is not a
covariate-adjusted spatial process model. Mixed spots may use this pattern test,
but spot variation also reflects mixtures and technical properties.

## Neighborhood enrichment and replicated comparisons

```shell
java -jar build/cli/jlinalg-0.3.6.jar spatial --method neighborhoods --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --units micrometer --coordinate-system synthetic-v1 --resolution cell --radius 10.1 --permutations 999 --seed 19 --out build/spatial-demo/neighborhoods
java -jar build/cli/jlinalg-0.3.6.jar spatial --method compare --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --units micrometer --coordinate-system synthetic-v1 --resolution cell --radius 10.1 --reference control --tested case --out build/spatial-demo/compare
```

Neighborhood tests permute cell-type labels within the same exchangeability
groups. This preserves observed label counts within each stratum and graph
geometry/boundaries. Each unordered type pair receives its observed edge count,
permutation mean, enrichment and depletion p-values. Both tails, all pairs and all
sections enter one BH family. Labels absent from a section remain in the family.
Zero-edge sections have no test. These tests require cell-resolved inputs; a
dominant cell-type label on a mixed spot is not a measured cell identity.

The replicated comparison first pools edge counts from all sections of each
biological sample. For every type pair it computes the fraction of all sample
edges belonging to that pair, transforms it by `asin(sqrt(fraction))`, and fits
the donor-aware sample-level model. Extra sections or edges never increase the
number of independent samples. Samples must have at least one edge. Use
`--paired true` with the paired sample sheet for complete donor pairs; numeric
covariates and batch handling follow the single-cell workflow.

This effect is relative edge composition. It can change because of population
abundance, tissue geometry or sampling, and is not an abundance-adjusted interaction
coefficient. Gaussian sample-level errors are assumed. Spatial autocorrelation is
kept inside the sample summary rather than supplying thousands of independent
edge-level observations to a condition test. Proximity and expression compatibility
do not establish causal signaling.

## Domains, niches and gradients

Pathology/domain labels can be supplied in a cell metadata column. Use
`single-cell --method state --group-column domain` to aggregate RNA by sample and
domain, or `single-cell --method abundance --group-column domain` for relative
domain abundance. `--method pathways --group-column domain` provides mean measured
gene-set expression. This is an annotated-domain workflow; exploratory domain
discovery, stability assessment and cross-sample domain registration are not yet
implemented. Labels must mean the same thing across samples. Testing the same
features used to discover a domain is exploratory without independent validation.

For a prespecified distance to an anatomical structure, provide a nonnegative
cell-level distance column in the declared physical units:

```shell
java -jar build/cli/jlinalg-0.3.6.jar spatial --method gradient --counts examples/cell-spatial/counts.tsv --cells examples/cell-spatial/cells.tsv --samples examples/cell-spatial/samples.tsv --features examples/cell-spatial/features.tsv --units micrometer --coordinate-system synthetic-v1 --resolution cell --distance-column distance --feature-list examples/cell-spatial/selected-features.tsv --reference control --tested case --out build/spatial-demo/gradient
```

The command calculates a descriptive linear slope of normalized expression versus
distance within each biological sample, pooling that sample's sections. It then
compares sample slopes using the donor-aware condition model; `sample-slopes.tsv`
contains the intermediate values. There are no spot-level significance claims.
Each sample needs >=3 retained observations and nonconstant distances. Unequal
precision of sample slopes is not modeled, distances/boundaries are treated as
known, and condition effects can reflect other distance-correlated cell mixtures.
Nonlinear trends, section-specific slopes and uncertainty in interfaces need a
different model and remain open.

Slope calculations center and scale both inputs. Constant responses have slope
zero; a slope outside the representable numeric range requires rescaling physical
units. Sample-level uncertainty follows the zero-residual-variance rule in the
[single-cell tutorial](single-cell.md).

## Limits and reproducibility

The shared input bounds apply. Graph construction is bounded at 25 million
within-compartment distance evaluations and one million edges. Permutation work
is bounded at 100 million observation/edge/pair visits per invocation; neighborhood
tests allow at most 50 cell labels. Gradient evaluation allows 20 million
observation-feature visits. Graph search is bounded quadratic, not an atlas-scale
spatial index. The counts stay sparse until a requested operation requires a
bounded representation or aggregate.
The one-million-edge cap also applies to supplied adjacency. Hypothesis/summary
tables allow 250,000 rows; donor-level comparisons follow the shared design and
total fitting-work limits. Radius search evaluates each unordered pair once;
kNN retains only its best k candidates, and section operations reuse partitioned
edges. Permutations reuse centered values and their fixed sum of squares.

Every run writes QC, original metadata, hashes, effective configuration, graph
parameters and null/estimand definitions. Seeds reproduce Java permutation results
for the same inputs. Record label definitions, physical scales, panel coverage,
upstream segmentation and registration versions, and sensitivity analyses.

Mixed-spot reference deconvolution, RNA/protein/ATAC integration, histology features,
image registration, 3D reconstruction, ligand-receptor activity inference and
uncertain-boundary modeling are still open [roadmap items](../../TODO.md).

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [P. A. P. Moran (1950) — Notes on continuous stochastic phenomena](../CITATIONS.md#moran-1950) — [PMID: 15420245](https://pubmed.ncbi.nlm.nih.gov/15420245/)
- [R. C. Geary (1954) — The Contiguity Ratio and Statistical Mapping](../CITATIONS.md#geary-1954)
- [Giovanni Palla et al. (2022) — Squidpy: a scalable framework for spatial omics analysis](../CITATIONS.md#palla-squidpy-2022) — [PMCID: PMC8828470](https://pmc.ncbi.nlm.nih.gov/articles/PMC8828470/)
- [Jordan W. Squair et al. (2021) — Confronting false discoveries in single-cell differential expression](../CITATIONS.md#squair-pseudoreplication-2021) — [PMID: 34584091](https://pubmed.ncbi.nlm.nih.gov/34584091/)
- [Yoav Benjamini and Yosef Hochberg (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](../CITATIONS.md#benjamini-hochberg-1995)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
