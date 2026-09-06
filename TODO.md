# Development inventory

Last reviewed: 2026-09-05.

This inventory keeps unfinished decisions at the top. Completed work is archived
at the bottom so it remains searchable without obscuring the active backlog.

Priority meanings:

- **Major**: removes a material scale limit or adds a cross-cutting analysis
  path.
- **Medium**: adds an important inference, workflow, or documentation feature
  without changing the core architecture.
- **Candidate**: a documented boundary that still needs a product decision.

## Active committed work

No major or medium items are currently in progress.

## Candidates requiring prioritization

These are real documented boundaries, but the vignettes do not establish that
they are committed work:

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
- **GLMM:** exact or adaptive-quadrature alternatives where PQL or first-order
  Laplace is inadequate.
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

## Completed work

### Frequentist zero-inflated mixed and pedigree models — completed

Completed 2026-09-05. `SparseZeroInflatedMixedModel` now fits ZIP and NB2-ZINB
models with sparse random effects in both the count and structural-zero
predictors. For zero observations the random-mode system includes the full
observed count, zero, and cross curvature. Independent terms accept identity or
caller-supplied precision matrices; `CorrelatedZeroInflatedRandomEffect`
estimates two variances and a guarded correlation with joint precision
`inverse(G) kron inverse(A)`.

BOBYQA remains derivative-free and sees only the low-dimensional outer
parameters. Sparse damped Newton solves the random mode. The optional
`fitWithInference` path computes a post-fit numerical observed Hessian of the
marginal Laplace objective, including nuisance covariance parameters, and
reports fixed/count/zero/dispersion covariance and standard errors. Prepared
fits reuse one symbolic sparse analysis and one numerical factor per worker;
profile-likelihood grids and deterministic conditional parametric bootstrap
reuse that prepared structure. Results expose boundary/singularity warnings.

Completion evidence:

- analytic two-predictor scores and observed Hessian blocks match finite
  differences for ZIP and ZINB zero and positive observations;
- ordinary and pedigree effects work independently in either process, and a
  simulation-recovery test resolves structural-zero group ordering;
- compiled count/zero/size formulas avoid repeated parsing and design builds;
- grouped estimates, marginal log likelihoods, variance components, and ZIP
  standard errors pass checked-in `glmmTMB` 1.1.14 fixtures;
- automatic bounded BFGS uses parallel full-Laplace numerical gradients for
  small modes while large sparse modes and failed BFGS fits retain BOBYQA;
- a separately compiled TMB 1.9.25 sparse-GMRF template gates correlated
  pedigree fixed effects, variances, correlation, and marginal likelihood;
- the optimized reproducible 500-member, 1,500-observation benchmark used
  1,000 random coefficients, 4,900 sparse equation/factor nonzeros, converged
  in 2.161532 seconds (previously 2.576353), and sampled a 312,168,592-byte
  peak-heap delta on the documented development host.

Sources: [pedigree and GLMM vignette](docs/vignettes/pedigree-and-glmm.md),
[`SparseZeroInflatedMixedModel`](src/main/java/org/jlinalg/distributional/SparseZeroInflatedMixedModel.java),
and the [R/TMB generators](src/test/resources/r-reference).

### Sparse Cox mixed and pedigree solves — completed

Completed 2026-09-05. Generic sparse precision, repeated observations,
prepared symbolic factorization, sparse diagnostics, dense/sparse parity, and a
deterministic `coxme` performance gate are implemented. The 800-row reference
scan was 10.03x faster with eight JLinAlg workers than single-threaded R.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Scalable pedigree prediction uncertainty — completed

Completed 2026-09-05. Sparse LMM diagonal PEV/reliability extraction uses
batched inverse-column solves, supports selected members, agrees with the dense
reference, and crosses the former 256-coefficient cutoff. Full dense covariance
remains intentionally quadratic rather than part of this contract.

Source: [pedigree vignette](docs/vignettes/pedigree-and-glmm.md).

### Prepared Cox score-scan pipeline — completed

Completed 2026-09-05. One null model and risk-set plan support deterministic
bounded, file-backed scans with model-based, cluster-robust, or supplied
correlation variance. The checked-in fixture measured 28,748 predictors/second.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### xWAS all-pairs output and BH/FDR correction — completed

The xWAS pipeline can stream all analyzable pairs, applies disk-backed BH over
the complete prespecified finite-p family, and records exclusions, failures,
test counts, and adjusted values explicitly.

Source: [xWAS pipeline](docs/vignettes/xwas-mr-pipeline.md).

### xWAS hit hand-off and robust follow-up — completed

Completed 2026-09-05. Retained hits preserve harmonized instruments and feed
RAPS, contamination-mixture, and PRESSO-style follow-up with deterministic
ordering, warnings, partial results, and CLI persistence.

Source: [xWAS follow-up](docs/vignettes/xwas-mr-pipeline.md#follow-up-analyses-for-hits).

### Dedicated colocalization workflow and MR documentation — completed

Completed 2026-09-05. Markdown and website vignettes cover alignment, priors,
H0-H4, shared-variant posterior, diagnostics, and MR/xWAS hand-off, with
navigation and stale statements repaired.

Source: [colocalization vignette](docs/vignettes/colocalization.md).

### Cox inference and diagnostics — completed

Completed 2026-09-05. Cluster sandwich covariance, residual/influence exports,
term/global proportional-hazards tests, and fixed/profiled gamma frailty are
available and independently gated.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Start-stop Cox performance path — completed

Completed 2026-09-05. An indexed ascending sweep maintains delayed-entry risk
moments while preserving Breslow/Efron behavior and the reference path.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Incremental omics result sink — completed

Confirmed 2026-09-05. Association and Cox scans stream deterministic estimates,
failures, and accounting rather than retaining large result grids in memory.

Source: [GWAS/TWAS pipeline](docs/gwas-twas-pipeline.md).

### Zero-inflated outer-point cache optimization — completed

Completed 2026-09-05. Worker-local caches remove repeated fixed-predictor,
NB2 gamma/size, and eligible zero-probability calculations from inner mode
solves. Full Salamanders fits improved another 1.41x (ZIP) and 1.83x (ZINB),
with identical reported estimates and standard errors. Added reproducible
Java/R runners, raw paired timings, independent likelihood reconstruction,
and parallel/cache-refresh tests. The large pedigree case showed no speedup.

Source: [benchmark and validation report](docs/zero-inflated-cache-performance.md).
