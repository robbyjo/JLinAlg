# Development inventory

Last reviewed: 2026-09-04.

This inventory collects explicit future-work statements from the documentation
and website vignettes, plus gaps found while reviewing the implemented APIs. It
tracks product priority separately from estimated effort: a small change can be
important, and a large change can remain deferred.

Priority meanings:

- **Major**: removes a material scale limit or adds a cross-cutting analysis
  path.
- **Medium**: adds an important inference, workflow, or documentation feature
  without changing the core architecture.
- **Candidate**: documented boundary that needs a product decision before it is
  promoted to committed work.

## Major priority

### Sparse Cox mixed and pedigree solves

The fixed right-censored Cox pass is already sorted and linear after setup, but
mixed-model and pedigree random information remains dense. Add a sparse
precision/factorization path for Gaussian shared, GRM, and pedigree frailty so
large related cohorts do not require a dense random-effect information matrix.

Done when:

- the sparse path does not materialize a dense random-effect covariance or
  information matrix;
- dense and sparse estimates agree on numerical fixtures;
- singularity, convergence, and backend behavior are reported explicitly; and
- a large-pedigree benchmark records time, peak memory, and numerical parity.

Sources: [Cox vignette](docs/vignettes/cox-survival.md),
[website survival vignette](site/vignettes/survival.html).

### Scalable pedigree prediction uncertainty

Provide prediction-error variance and reliability from the sparse combined
pedigree model. The API and documentation should use population-neutral terms
such as individuals or pedigree members; animal breeding is one use case, not
the only one.

A complete conditional covariance is inherently quadratic in the number of
predicted individuals. The scalable target should therefore be diagonal PEV and
reliability for all requested individuals, with optional selected covariance
blocks, rather than promising that an arbitrarily large full dense matrix is
cheap.

Done when:

- diagonal PEV/reliability can be extracted without forming the full inverse;
- callers can request selected individuals or covariance blocks;
- results agree with the dense reference path on small pedigrees; and
- memory scaling is documented and benchmarked.

Sources: [pedigree vignette](docs/vignettes/pedigree-and-glmm.md),
[website pedigree vignette](site/vignettes/pedigree.html), and the
[compatibility roadmap](docs/lme4-pedigreemm-roadmap.md).

### Prepared Cox score-scan pipeline

Add a prepared null-model score scan for many genetic or molecular predictors,
including file-backed blocks and the existing cohort/QC conventions. This is
the survival-analysis analogue of the association pipelines and should reuse a
fixed Cox risk-set plan rather than refitting the null model for every marker.

Done when:

- right-censored score tests reuse one fitted null model and risk-set plan;
- block streaming, deterministic output, failures, and provenance match the
  association pipeline conventions;
- robust and relatedness-aware variance can be selected when those engines are
  available; and
- numerical fixtures and a representative high-throughput benchmark are
  checked in.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

## Medium priority

### xWAS all-pairs output and BH/FDR correction — completed

The xWAS MR pipeline retains threshold-passing pairs for full diagnostics and
can now stream every successfully screened pair to an optional all-pairs table.
The CLI appends a BH-adjusted value across the complete analyzable family.

The implementation reuses the CLI's disk-backed `ExternalBh` path, so the
bounded scan does not retain the complete exposure-by-outcome grid in memory.

Completed behavior:

- BH uses every finite screening p-value in the prespecified exposure-by-outcome
  family, without a nominal-p prefilter;
- excluded, insufficient-instrument, and failed pairs have explicit semantics
  and are not silently counted as tests;
- the CLI records the number of adjusted tests and emits `fdr_bh`; and
- adjusted values match an independent R/JDistlib reference including ties,
  zero, one, and non-finite results.

Sources: [xWAS pipeline vignette](docs/vignettes/xwas-mr-pipeline.md) and
[`ExternalBh`](src/main/java/org/jlinalg/cli/ExternalBh.java).

### xWAS hit hand-off and robust follow-up

Make retained hits easy to pass into RAPS, contamination-mixture, and
PRESSO-style analyses without redoing the full scan. Prefer a bounded design:
retain or reconstruct harmonized instruments only for hits rather than storing
them for every pair. Decide whether follow-up runs automatically or through a
second explicit CLI stage, and persist convergence, outliers, warnings, and
method-specific p-values.

Source: [xWAS hit follow-up](docs/vignettes/xwas-mr-pipeline.md#follow-up-analyses-for-hits).

### Dedicated colocalization vignette and MR documentation repair

`ColocSusie` is implemented, tested against `coloc::coloc.susie`, benchmarked,
and briefly demonstrated inside the combined SuSiE/SEM vignette. What is
missing is a dedicated colocalization workflow. Add a vignette covering input
alignment, priors, posterior overlap/trimming, H0-H4 interpretation,
variant-level shared posterior, diagnostics, and an MR/xWAS follow-up example.

At the same time, remove the stale end-to-end MR statement that JLinAlg does not
implement colocalization, add the new page to both vignette indexes/navigation,
and keep the generated website synchronized with the Markdown source.

Sources: [current compact example](docs/vignettes/susie-and-sem.md),
[stale MR statement](docs/vignettes/mr-end-to-end.md), and
[colocalization implementation](src/main/java/org/jlinalg/coloc/ColocSusie.java).

### Cox inference and diagnostics

Add the following as separable deliverables so they can land independently:

1. cluster-robust sandwich covariance, including recurrent-event use;
2. proportional-hazards diagnostics based on Schoenfeld/score residuals and a
   documented global/test-per-term interface;
3. influence and residual exports suitable for diagnostic plots; and
4. gamma frailty, with the approximation and integration method explicit.

These should be validated independently rather than hidden behind one broad
"Cox diagnostics" completion flag.

Sources: [Cox vignette](docs/vignettes/cox-survival.md),
[website survival boundary](site/vignettes/survival.html).

### Start-stop Cox performance path

Start-stop observations currently use the general risk-set reference path.
Profile representative recurrent-event and time-dependent-covariate workloads,
then add an indexed/sweep implementation if the benchmark confirms a material
bottleneck. Preserve the reference path for parity tests.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Incremental omics result sink

Complete the planned incremental sink for omics response scans so large
feature-by-predictor results need not be retained in memory. Match the existing
association sink's accounting, deterministic ordering, and failure contract.

Source: [GWAS/TWAS pipeline guide](docs/gwas-twas-pipeline.md).

## Candidates requiring prioritization

These are real documented boundaries, but the vignettes do not establish that
they are committed major- or medium-priority work:

- **SEM:** latent measurement variables, mean/intercept structures, ordinal
  thresholds, robust/clustered corrections, modification indices,
  indirect-effect delta-method inference, and FIML missing-data patterns.
- **Meta-analysis:** cluster-robust variance, multilevel/correlated effects,
  publication-bias diagnostics, and effect-size construction.
- **Mixed models:** sparse correlated-block likelihoods, scalable sparse
  finite-DF calculations, broader `lme4` formula parity, and profile-likelihood
  intervals.
- **Time series:** diffuse exact likelihood for integrated models and a sparse
  missing-series path.
- **MR workflow:** a general estimator CLI, native plot rendering, and
  conditional/secondary-signal clumping.
- **GLMM:** exact or adaptive-quadrature alternatives where PQL is inadequate.
- **Penalized regression:** selection-aware inference after LASSO/elastic-net
  selection.

Sources: [SEM](docs/vignettes/sem.md),
[meta-analysis](docs/vignettes/meta-analysis.md),
[mixed-model roadmap](docs/lme4-pedigreemm-roadmap.md),
[time series](docs/vignettes/time-series.md),
[end-to-end MR](docs/vignettes/mr-end-to-end.md),
[GLMM/pedigree](docs/vignettes/pedigree-and-glmm.md), and
[linear/penalized models](docs/vignettes/linear-models-and-glms.md).

## Boundaries that are not automatically TODO items

Some vignette warnings state a statistical contract rather than a missing
feature. Do not turn these into implementation work without a separate design
decision. Examples include PQL not being an exact marginal GLMM likelihood,
ordinary post-selection OLS p-values not being selection-adjusted, MR
sensitivity analyses not proving the exclusion restriction, and the xWAS
benchmark not measuring an end-to-end TwoSampleMR workflow.
