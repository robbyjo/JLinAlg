# Single-cell and spatial workflow validation

The source build provides `single-cell` and `spatial` workflows for **quantified,
annotated RNA counts**. These are a tested first delivery of the large roadmap,
not completion of all protein, epigenetic, trajectory or imaging analysis suites.
The [single-cell tutorial](vignettes/single-cell.md) and
[spatial tutorial](vignettes/spatial-analysis.md) specify schemas, estimands and
unsupported designs. Unimplemented branches remain explicit in [TODO.md](../TODO.md).

## Delivered operations

| Operation | Scientific unit and supported behavior |
|---|---|
| QC | Sparse raw RNA counts, library/detection/mitochondrial metrics, imported exclusions, complete metadata and measured panel |
| Pseudobulk | Exact counts summed by biological sample and supplied population/domain; cells and sections remain nested |
| Differential state | R limma voom, independent donors or complete two-condition pairs, numeric covariates and batch, joint feature/population BH |
| Relative abundance | Equal-sample Gaussian inference for arcsine-square-root population proportions with explicit denominators |
| Functional scores | Prespecified measured gene-set mean log2 CPM and sample-level contrasts |
| Representation | Library normalization, variance selection, centered PCA and seeded k-means; exploratory |
| Spatial graph | Physical 2D/3D radius, gap-capped union kNN or supplied adjacency; section/compartment restrictions; SVG maps |
| Autocorrelation | Within-section Moran I / Geary C and restricted random-label permutations |
| Neighborhoods | Label-pair enrichment/depletion with anatomical strata and fixed geometry |
| Replicated relationships | Edge-pair fraction per sample; independent or paired donor comparison |
| Spatial gradients | Descriptive sample-specific linear slopes followed by donor-level condition comparison |

All commands preserve supplied annotations and hierarchy, write hashes and effective
configuration, and stage a new output directory before publishing. No implicit
package installation, remote annotation, estimator fallback or embedded LLM is used.

## Reproduce the independent comparisons

```shell
python tools/generate_cell_spatial_examples.py
./gradlew test --tests '*SpatialStatisticsTest' --tests '*SampleInferenceTest' --tests '*CellSpatialCliTest' executableJar
python tools/validate_cell_spatial.py --rscript Rscript --r-library /path/to/R/library
./gradlew check executableJar
```

On Windows use `gradlew.bat`; Python requires NumPy, and R requires limma. The
validation script creates a fresh timestamped `build/cell-spatial-*` directory,
with every command's console log, outputs, runtime, independent reference script,
reference results and R session. It never overwrites an existing audit directory.

The initial 2026-09-20 validation used Java on Windows 11, R 4.6.1, limma 3.68.5
and statmod 1.5.2. Fifteen end-to-end invocations passed, covering independent and
paired state/abundance, all three graph modes, and every documented operation.
The fixture has 320 observations, 8 samples, 16 sections, 31 measured genes and
9,920 nonzero entries. Each invocation took approximately 0.22–0.82 seconds on the
tested host, including process startup; this is a small example, not an atlas-scale
throughput or memory benchmark.

The repository-wide gate discovered **919 tests: 916 passed, 0 failures/errors,
3 optional CHOLMOD skips**. The executable JAR, Javadoc, website structural/link
checks and generated-citation consistency checks also passed. The ordinary full
test run stalled in `OpenClComputeBackend.initialize` / `CL.clBuildProgram`, as
confirmed by a JVM thread dump. It was stopped and the suite rerun with a temporary
test-classpath copy of JDistlib omitting only
`META-INF/services/jdistlib.accelerator.ComputeBackend`, so optional accelerator
discovery could not block CPU validation. Production dependency bytes and the
executable JAR retain their normal accelerator configuration. This run therefore
does not claim validation of those optional GPU backends.

Independent checks reconstruct the count aggregation directly from sparse raw
inputs. They refit both independent and paired state models directly in R from
raw counts and metadata, using `xtabs` and `model.matrix` rather than the Java
aggregate/design files. Effects, SEs and p-values agree within `rtol=1e-8` and
`atol=1e-11`. Agreement checks adapter correctness and design alignment, not
biological validity of the model assumptions.

Abundance coefficients and SEs also match a separate NumPy least-squares solution;
R `lm` supplies an independent check of Student t p-values. Moran/Geary values
match an independent dense adjacency-matrix calculation for every requested
section/feature (`rtol=1e-10`, `atol=1e-12`). Monte Carlo permutation results are
seed reproducible. The production command and reference both call limma, so this
is not an independent reimplementation of limma's estimator.

## Statistical and failure fixtures

- The four-node path has independently calculated Moran I = 1/3 and Geary C = 0.3.
  All 24 assignments are enumerated separately; seeded Monte Carlo tails agree
  with the exhaustive distribution within sampling tolerance.
- Independent two-group and complete-pair estimates and uncertainty match hand
  calculations. A 400-replicate iid Gaussian null checks broad calibration of the
  reused sample-level model; it does not establish calibration for rare cell types,
  heterogeneous sequencing depth or biological single-cell datasets.
- Duplicate observations/features/count entries, invalid counts, missing IDs,
  unknown options and inappropriate per-method arguments reject.
- Repeated donors in independent designs, incomplete pairs, condition/batch
  confounding, insufficient residual degrees of freedom and mixed sources/assays
  reject or receive an explicit unsupported-population status.
- Count totals survive aggregation. Zero-library and high-mitochondrial cells
  retain exclusion reasons. The full assayed feature universe survives QC.
- Overlapping numeric coordinates never connect independent sections. Three-
  dimensional distance, gap limits, compartment restrictions, isolates, constant
  features and singleton permutation strata are exercised.
- Mixed-spot cell-label neighborhood inference rejects. Existing outputs are
  preserved; failed jobs do not publish a partial result directory.

## Correctness and performance audit, 2026-09-20

The follow-up audit reproduced and fixed these defects against the first delivery:

- R type inference changed numeric-looking sample/observation IDs and literal
  `NA` feature IDs; unquoted output could damage embedded quotes. Explicit character
  parsing and quoted UTF-8 output now preserve identities. The old adapter failed
  its count/design alignment check on leading-zero sample IDs. Java also verifies
  the returned feature schema and ordering.
- An exactly fitted native sample response produced a spurious p-value near
  `1.4e-91` from numerical residuals. Residual RMS at the floating-point floor now
  produces unavailable inference. Constant decimal spatial values, extreme scales
  (`1e-200`, `1e200`, and opposite-sign `1e308`), and tiny physical units have
  explicit regression coverage; centered/scaled products preserve the estimand.
- Chained distance rounding omitted some three-dimensional pairs exactly on an
  inclusive radius boundary. Independent all-pair fixtures now cover radius and
  union kNN, ties, compartments and small distances.
- Partially donor-nested batch indicators could reject identifiable paired
  designs. Batch basis selection now removes only dependence on donor effects and
  other batch indicators; genuine condition/batch confounding still rejects.
- Invalid covariates and gene sets are validated before population eligibility
  can mask an error. Empty test families receive accurate status and group-directory
  mappings. Dense designs, output tables, total fitting work and PCA now have
  explicit limits before costly allocations or decompositions.

The audited version passed **18 end-to-end invocations**, including identifier
round trips, unchanged numerical results after renaming, and all-filtered genes.
The repository-wide gate found **929 tests: 926 passed, no failures/errors, and
3 optional CHOLMOD skips**. Javadoc and executable packaging passed. The full
suite uses the same CPU validation workaround described above; optional GPU
backends remain outside this validation.

Reproduce the standalone kernel benchmark with:

```shell
python tools/benchmark_cell_spatial.py --javac /path/to/javac
```

The script compiles frozen baseline commit `3a8b4f12eee31db935efbada15150bdbd94e9830`
and current spatial code in separate JVMs, checks graph counts/checksums and numerical
agreement, and saves raw timings and medians. The following host measurements use
Java 25, three warmups, seven measurements and `-Xms256m -Xmx1g`:

| Workload | Baseline median | Revised median | Ratio |
|---|---:|---:|---:|
| Radius: 600 nodes, 179,700 edges | 183.10 ms | 7.49 ms | 24.45x |
| Union kNN: 600 nodes, k=6, dense candidates | 36.53 ms | 5.21 ms | 7.01x |
| 5,000 permutations: 800 nodes, 799 edges | 41.01 ms | 30.49 ms | 1.34x |

These are bounded synthetic kernel timings, not general end-to-end speedups or
atlas-scale memory/throughput validation. Timing is reported rather than asserted
as a test threshold and varies with host load and JVM behavior. Section edge
partitioning, cached feature lookups and gene-set
parsing also remove repeated scans without changing scientific units.

## Open scientific and engineering work

The tests do **not** constitute a public multi-donor acceptance study. Real-cohort
replication, RNA null/FDR/coverage under heterogeneous donors, rare-population
stress tests, imperfect references, extensive scale/memory benchmarks and
annotation/segmentation sensitivity remain necessary. Current constraints and
outcomes are deliberately distinguished from these unperformed validations.

The broader roadmap still includes automatic ambient/doublet assessment, reference
annotation and marker discovery, neighborhood DA/Milo, multi-phenotype covariance,
CITE-seq/totalVI, mass-spectrometry/QFeatures, cytometry, ATAC/Signac, trajectories,
velocity, perturbations, regulatory/communication models and integrated genetic
follow-up. Spatial remainder includes deconvolution/cell2location, domain discovery,
nonlinear covariate-adjusted spatial processes, image features, registration,
uncertain interfaces, multi-section reconstruction and multimodal reporting.

These are not implemented by merely exporting tables or naming a candidate adapter.
The native table paths are not claimed equivalent to complete Scanpy, Squidpy,
muscat or Seurat installations. Library-size normalization is not TMM; pathway
mean-expression scores are not regulator activity; spatial edge composition is
not causal signaling; extra cells/sections are not extra independent donors.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Helena L. Crowell et al. (2020) — muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data](CITATIONS.md#crowell-muscat-2020) — [PMID: 33257685](https://pubmed.ncbi.nlm.nih.gov/33257685/)
- [Jordan W. Squair et al. (2021) — Confronting false discoveries in single-cell differential expression](CITATIONS.md#squair-pseudoreplication-2021) — [PMID: 34584091](https://pubmed.ncbi.nlm.nih.gov/34584091/)
- [Gordon K. Smyth (2004) — Linear models and empirical Bayes methods for assessing differential expression in microarray experiments](CITATIONS.md#smyth-limma-2004) — [PMID: 16646809](https://pubmed.ncbi.nlm.nih.gov/16646809/)
- [Charity W. Law, Yunshun Chen, Wei Shi, and Gordon K. Smyth (2014) — voom: precision weights unlock linear model analysis tools for RNA-seq read counts](CITATIONS.md#law-voom-2014) — [PMID: 24485249](https://pubmed.ncbi.nlm.nih.gov/24485249/) · [PMCID: PMC4053721](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053721/)
- [Yoav Benjamini and Yosef Hochberg (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](CITATIONS.md#benjamini-hochberg-1995)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
