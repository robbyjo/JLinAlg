# Unreleased

## Statistical audit repairs — 2026-09-08

- Added joint latent/mean/ordinal SEM and direct constrained FIML with full
  covariance, robust inference, modification indices, and indirect effects.
- Added multilevel random-moderator ML/REML covariance and correlated CR2
  inference; replaced median reflection with L0/R0/Q0 trim-and-fill.
- Extended sparse unstructured mixed models, finite-DF calculations, formula
  handling, and nuisance-optimized profile refits.
- Replaced diffuse bookkeeping and dense missing-series covariance with
  exact diffuse, state-sized Kalman filtering and filtered-state forecasts.
- Separated conditional instrument associations from causal MR inference;
  added conditional signal selection, CLI controls, and corrected SVG fits.
- Replaced fixed-node GLMM integration with mode-adaptive, node-refined
  quadrature, boundary checks, trials/offsets, and covariance inference.
- Fixed quantile gradients, multinomial convergence, kernel underflow, and
  supersmoother edge behavior; reused multivariate QR, added smoother-aware
  Robinson covariance, and added a certified nonsmoothed quantile LP.
- Preserved weights in held-out selection inference and added fixed-penalty
  Gaussian polyhedral LASSO/elastic-net inference.
- Added independent R fixtures and paired accuracy/speed measurements, with
  explicit slower cases and remaining limits. See
  [the validation report](docs/advanced-validation.md).

## Earlier additions

- Added sparse first-order Laplace beta mixed models with analytic beta
  Fisher-scoring updates, jointly optimized constant precision and variance
  components, grouped random effects, arbitrary sparse coefficient precision,
  and a pedigree convenience API consuming sparse `A^-1`. An independent
  `glmmTMB` 1.1.14 fixture checks fixed effects, precision, variance, likelihood,
  unrelated-founder equivalence, and retention of unphenotyped ancestors. On
  the documented 30,000-row/1,000-group workload, the Java median was
  1,028.105 ms versus 11,520 ms for R (11.20x faster); a 3,000-member pedigree
  fit took 1,374.480 ms.
- Added specialized classical beta regression compatible with R `betareg`'s
  mean/precision parameterization, six mean links, three precision links,
  constant and variable precision, expected-information covariance, and Wald
  inference. Package-example fixtures cover GasolineYield and FoodExpenditure;
  a deterministic benchmark compares the fused row-major solver with the
  generic distributional engine. On the documented 100,000-row workload,
  JLinAlg matched all coefficients within `1.3e-14` of its generic engine and
  within `2.6e-13` of R, while running in 143.216 ms versus 2,090 ms for
  `betareg` 3.2-6 (14.59x faster).
- Further reduced ZIP/ZINB likelihood work with worker-local outer-point
  caches for fixed predictors, NB2 size/gamma terms, and eligible zero-process
  probabilities. Fresh paired Salamanders full-fit medians improved from
  32.290 to 22.861 ms (1.41x) for ZIP and 211.470 to 115.625 ms (1.83x) for
  ZINB, with identical reported estimates/SEs and unchanged optimizer controls.
  Added reproducible Java/R benchmarks and cache/parallel correctness tests;
  the large correlated-pedigree benchmark showed no speedup. See
  [measurement details](docs/zero-inflated-cache-performance.md).
- Reduced sparse ZIP/ZINB fit time with reusable likelihood workspaces,
  response-only caches, likelihood-only Newton line searches, split
  count/zero sparse products, and bounded BFGS over workload-gated parallel
  full-Laplace numerical gradients. Automatic selection retains BOBYQA for
  large sparse modes and as a convergence fallback. On the official 644-row
  `glmmTMB` Salamanders examples, median JLinAlg fit time fell from 0.135 to
  0.033 seconds for ZIP and from 2.365 to 0.212 seconds for ZINB; both fits
  converged to the documented `glmmTMB` likelihoods.
- Added a uniform `org.jlinalg.stats.StatisticalTests` facade for the basic
  base R test families. Existing JDistlib tests are delegated; JLinAlg adds
  correlation, Fisher exact, Friedman, Cochran-Mantel-Haenszel, McNemar,
  one-way ANOVA, ANOVA power, proportion/trend, and Quade implementations.
- Added deterministic all-screened-pairs xWAS MR output with disk-backed
  Benjamini-Hochberg adjustment through `mr-xwas --fdr-output`.
- Corrected MR and SuSiE documentation to reflect the implemented
  `ColocSusie` workflow.

# JLinAlg 0.2.0

JLinAlg 0.2.0 adds independently validated, performance-tuned structural
equation modeling, complete end-to-end Mendelian-randomization workflows, and
SuSiE colocalization and reference parity.

## Structural equation models

- Replaced derivative-free SEM optimization and finite-difference inference
  with an analytic RAM likelihood gradient, scaled BFGS optimization, and
  analytic expected Fisher information.
- Corrected normal-theory ML chi-square/RMSEA sample-size scaling and aligned
  Bentler SRMR standardization with `lavaan`.
- Added a `lavaan` 0.7-2 regression fixture covering paths, variances,
  exogenous and residual covariances, standard errors, likelihood, and fit
  indices.
- Added a 4,680-row, 12-variable TOPMed cardiometabolic benchmark. JLinAlg and
  `lavaan` agreed below `1e-8` for parameter estimates and below `5e-10` for
  standard errors; warmed median fit times were 0.00816 and 0.0700 seconds,
  respectively (8.58x faster on the documented host).
- Added a standalone SEM vignette and public verification/results webpage.

## Mendelian randomization and LD

- Added end-to-end instrument discovery, download, formatting, PLINK LD
  clumping, and xWAS MR CLI workflows.
- Added bounded parallel exposure-by-phenotype screening with reusable clumped
  instruments, scale-safe significance filtering, diagnostics for retained
  hits, explicit failure output, and reproducible metadata.
- Added a 45,000-pair Java/base-R benchmark with floating-point-roundoff
  agreement and a measured 5.92x eight-thread speedup.

## SuSiE and colocalization

- Added multi-signal `ColocSusie` with ID alignment, posterior-overlap
  trimming, scalar or weighted priors, H0-H4 summaries, and per-variant H4
  posteriors.
- Matched susieR sample scaling, finite-sample z transformation, and ELBO
  convergence using the official `N3finemapping` data.
- Removed repeated hot-loop allocations, parallelized symmetric
  cross-products, and cached per-effect `X'X b` updates.
- On the documented 574-by-1,001 benchmark, JLinAlg matched fixed-prior susieR
  results within `2e-10` and ran 4.83x faster.

## Verification and compatibility

- All shipped feature families remain covered by automated correctness tests
  and executable performance checks. Cross-language gates use base R,
  statsmodels, `nlme`, `rrBLUP`, `lme4`, `mgcv`, `gamlss`, `VGAM`,
  `MendelianRandomization`, `susieR`, `coloc`, and `lavaan` where comparable.
- The clean release gate discovered 253 tests: 250 passed, zero failed, and
  three optional CHOLMOD-native tests were skipped in the portable build.
- Published speed ratios are workload-specific measurements, not universal
  guarantees; the verification webpage links each number to its statistical
  contract and reproduction commands.
- Requires Java 17 or newer and JDistlib 0.10.1.

## Build

Run `./gradlew check executableJar` (or
`.\gradlew.bat check executableJar` on Windows). The self-contained artifact
is `build/cli/jlinalg-0.2.0.jar`.
