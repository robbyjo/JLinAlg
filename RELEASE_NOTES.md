# Unreleased

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
