# Development inventory

Last reviewed: 2026-09-08.

The September 8 audit found that several earlier completion claims described
approximations or scaffolding rather than the requested methods. Those claims
are superseded by this inventory and the
[validation report](docs/advanced-validation.md).

## Audited implementation status

The requested analysis paths and audited repairs are implemented and independently
rechecked. The final integrated gate passes 459 tests with three optional native
CHOLMOD skips. The remaining extensions below are explicit limits, not claims
that the requested methods have full R-package parity.

## Remaining extensions and limits

These boundaries are explicit; this release does not claim complete parity
with every feature of lavaan, metafor, lme4, or the other R packages.

- SEM: DWLS/WLSMV, mixed/missing ordinal responses, multigroup invariance,
  ordinal modification indices, and FIML robust scaled fit statistics.
- Meta-analysis: sparse large-study covariance estimation, variance-component
  profiles, and additional structured covariance families.
- Mixed models: joint multi-DF KR tests, nonregular boundary-profile coverage,
  general correlation profiles in blocks larger than two, original-coordinate
  full PEV blocks, and joint unstructured-plus-pedigree covariance optimization.
  Formula-level pedigree mapping and automatic complete-case alignment remain open.
- Time series: diffuse coefficient covariance, general regression terms, and
  historical state smoothing.
- MR: additional generalized estimators requiring richer summary statistics;
  conditional GWAS calculations currently use the documented Gaussian score model.
- GLMM: multidimensional/crossed/pedigree adaptive quadrature beyond the new
  scalar Gaussian random-intercept path.
- Selective inference: unknown-noise or response-selected penalties, dependent
  sample splits, and non-Gaussian selection require different inferential methods.
- Quantile/nonparametric regression: quantile inferential covariance,
  mixed-type/multidimensional kernels, and automatic inference-valid bandwidth
  selection remain distinct extensions; exact nonsmoothed quantile fitting is implemented.

## Completed audit repairs — 2026-09-08

- Joint latent RAM, structural means/intercepts, ordinal probit pairwise ML,
  direct constrained FIML, robust/cluster covariance, efficient continuous-model
  modification indices, and indirect-effect delta inference.
- Joint ML/REML multilevel random-moderator covariance, correlated CR0/CR1/CR2
  inference, PET/PEESE, real L0/R0/Q0 trim-and-fill, and stable effect sizes.
- Sparse multi-block unstructured covariance fits, model-derived finite DF,
  weights/offsets, and formula-native nuisance-optimized ML profiles.
- Exact diffuse integrated/seasonal likelihood and state-sized missing-series
  filtering and forecasting.
- Generalized/overlap-aware estimator CLI, native SVG fits, explicit conditional
  instrument p-values, and conditional forward/backward secondary-signal selection.
- Mode-adaptive, node-refined binomial/Poisson random-intercept quadrature,
  offsets/trials, zero-variance comparison, and covariance inference.
- Correct weighted held-out LASSO/elastic-net inference and exact Gaussian
  fixed-penalty active-set/sign polyhedral inference.
- Shared-QR multivariate OLS, score-checked multinomial optimization, corrected
  smoothed quantile gradient and certified nonsmoothed quantile LP, stable
  Gaussian kernel prediction, R-compatible supersmoother, and smoother-aware
  Robinson partially linear slope inference.
- Independent R fixture generators, regression tests, accuracy-gated benchmark
  runners, paired measured timings, source vignettes, and website documentation.

Sources: [advanced extensions](docs/vignettes/advanced-extensions.md),
[regression families](docs/vignettes/regression-families.md),
[mixed-model roadmap](docs/lme4-pedigreemm-roadmap.md), and
[accuracy and performance](docs/advanced-validation.md).

## Earlier completed work

### Nonlinear fixed and mixed effects

`NonlinearFixedModel` uses analytic gradients and damped Gauss–Newton.
`NonlinearMixedModel` linearizes additive Gaussian nonlinear means while
reusing sparse REML equations; `fitPedigree` adds numerator-relationship
precision terms. This is not an exact general `nlme` likelihood or a
non-Gaussian nonlinear quadrature model.
See [the vignette](docs/vignettes/nonlinear-models.md).

### Frequentist linear mediation analysis — completed

Completed 2026-09-07. `MediationAnalysis` fits the mediator, outcome, and
total-effect Gaussian OLS models and reports the `a`, `b`, indirect, direct,
and total effects. Sparse REML mediation also supports ordinary random-effect
terms and additive pedigree precision terms, reusing prepared sparse
structure across all three fits. Indirect-effect inference uses the analytic
Sobel/product-of-coefficients delta method with a normal confidence interval;
no bootstrap, Monte Carlo, or posterior sampling is used. Common complete-row
selection keeps all three component fits on the same analysis sample.

Source: [mediation vignette](vignettes/mediation.md).

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
