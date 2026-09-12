# Unreleased

## Meta-analysis CLI and array workflows

- Added `meta-analysis` and `meta-regression` commands for separate cohort
  TSV/CSV files (including gzip), matched by feature ID with external sorting.
- Fixed/random effects support REML, DL, and PM, numeric moderators, coefficient
  inference controls, cohort direction strings, and `--min-cohorts` (default 1).
  Single-cohort pooling passes through explicitly; insufficient or rank-deficient
  regression rows are marked as unestimable.
- Extended primitive-array batches to support missing cohorts, directions, and
  minimum counts; added primitive-array meta-regression. Heterogeneity stopping
  rules respect effect units and fail explicitly on iteration exhaustion.
- Added independent R/metafor CLI fixtures, runnable six-cohort examples, and a
  reproducible synthetic allocation/time comparison with the study-list API.

## Statistical correctness fixes

- Stabilized the Lugannani–Rice kernel-tail fallback at and around the mean,
  including quantile inversion. It remains an explicitly labeled approximation;
  invalid numerical probabilities now fail instead of being clamped.
- Made sparse REML and Paule–Mandel heterogeneity optimization independent of
  effect units and refined REML convergence near a flat numerical optimum.
- Enforced the initial multidimensional quadrature node budget before allocating
  dense random-effect matrices or fitting, including node-count overflow.
- Added scale-aware shape, symmetry, diagonal, and positive-semidefinite checks
  for conditional MR exposure covariances. Valid singular covariance matrices
  remain supported when conditional residual variances are positive.
- Applied the ordinary profile-likelihood input checks to boundary profiles,
  with boundary confidence levels restricted to `(0.5, 1)`.

Regression validation: 667 tests, 664 passed, zero failures, and three optional
native CHOLMOD tests skipped; website and Javadoc checks passed. Analytic
equal-variance meta-analysis and independent R integration checks were also run.

# JLinAlg 0.3.5

## Executable mediation and fine-mapping workflows

- The executable JAR now exposes Gaussian `mediation` workflows for ordinary,
  grouped-REML, and pedigree-REML models, with common complete-case alignment
  and explicit component-model convergence output.
- New `susie` and `coloc` commands run summary-statistic fine mapping from
  labeled LD matrices and multi-signal colocalization from the resulting
  effect tables.
- Progressive CLI tutorials now cover association, Mendelian randomization,
  mediation, SuSiE, and colocalization from input files through interpretation.

## Association preflight and option validation

- `--max-maf` and `--max-mac` complement the existing lower-bound variant
  filters and are applied after sample alignment.
- Phenotype-only `--dry-run` and `--explain` now compile and validate the
  selected route before fitting. Genotype transforms, unsupported family/link
  combinations, and the common `--explains` typo fail with specific guidance.

## Pedigree identity and hierarchy

- Family columns now disambiguate duplicate member IDs without partitioning
  the ancestry graph. Globally unique raw observation and parent IDs resolve
  automatically across family labels, while exact `family:individual` IDs
  remain accepted and ambiguous raw IDs fail with a qualification request.
- Parent-child links recursively define the full pedigree hierarchy, including
  shared ancestors and consanguineous matings. A direct R `pedigreemm` 0.3.5
  fixture now verifies the additive relationship matrix, inbreeding
  coefficients, and Henderson sparse inverse through a cousin-mating example.
- Runs with zero aligned source-pedigree matches now fail rather than silently
  constructing an all-singleton model. Logs and manifests distinguish matched
  source observations, automatically resolved aliases, and true singletons.

Release downloads: `jlinalg-0.3.5.jar` is the self-contained executable;
`JLinAlg-0.3.5-library.jar` is the thin library. Sources and Javadoc JARs and
`SHA256SUMS.txt` are also included.

# JLinAlg 0.3.4

## High-core adaptive omics scheduling

- Automatic feature blocks are sized from both current JVM heap headroom and
  the requested thread count.
- When memory permits, each block queues at least two complete 256-feature
  chunks per scan worker, giving work stealing enough queued work to absorb
  uneven feature-fit times without leaving cores idle at block barriers.
- Under tighter memory, JLinAlg selects the largest complete worker wave that
  fits or safely caps worker capacity. An explicit `--block-size` remains
  authoritative.

## Pedigree singleton handling

- Phenotype individuals absent from the pedigree file are retained as
  unrelated, noninbred singleton founders instead of stopping the analysis.
- Repeated observations with the same missing pedigree ID remain grouped in
  one singleton family. The console, log, and manifest report singleton
  families and observations.

## CLI documentation and transform semantics

- Added a comprehensive CLI-only guide and copy-pasteable command-line examples
  to the model vignettes, including grouped REML, pedigrees, GRMs, GLMMs,
  genotype scans, Cox models, and specialized subcommands.
- Omics row transforms intentionally run after phenotype complete-case
  filtering, so `winsor_mad` and other row-statistic transforms use the final
  analysis sample set.

Release downloads: `jlinalg-0.3.4.jar` is the self-contained executable;
`JLinAlg-0.3.4-library.jar` is the thin library. Sources and Javadoc JARs and
`SHA256SUMS.txt` are also included.

# JLinAlg 0.3.3

## Complete-case sample handling

- Phenotype rows with missing values in any model column are omitted by
  default, and the corresponding samples are removed from streamed omics
  inputs before fitting.
- Console, log, and manifest output distinguish aligned samples from the final
  analysis sample count and report how many incomplete phenotype rows were
  omitted.

## Exact streamed mixed-model scans

- Numeric Gaussian omics scans now refit variance components by REML for every
  feature through the `lmer`-like pathway instead of using P3D or EMMAX.
- Non-Gaussian numeric omics scans now refit a first-order Laplace GLMM for
  every feature through the `glmer`-like pathway instead of PQL.
- Pedigree-enabled scans use direct sparse Henderson relationship precision
  with ancestry-graph inbreeding, matching the intended `pedigreemm` model
  structure without materializing a dense relationship matrix.
- The selected mixed-fit strategy is recorded in the console preflight, run
  log, and manifest. Genotype LMM scans retain their documented null-model
  P3D pathway.

## Delimited output and compact schemas

- A `.csv` output path now produces correctly quoted comma-separated output;
  `.tsv` and other existing outputs remain tab-separated.
- Generic numeric omics results omit genotype-only allele, frequency, call,
  imputation, Hardy-Weinberg, and filter columns.
- Per-row error fields are consolidated as `failure_reason`. Invariant
  `omics_type`, `statistic_type`, `df_method`, `partial_r2_method`, and
  `output_format` metadata now live in the run log and manifest instead of
  being repeated in every result row.

Release downloads: `jlinalg-0.3.3.jar` is the self-contained executable;
`JLinAlg-0.3.3-library.jar` is the thin library. Sources and Javadoc JARs and
`SHA256SUMS.txt` are also included.

# JLinAlg 0.3.2

## Automatic phenotype/omics sample alignment

- Omics scans now intersect phenotype and omics sample IDs automatically,
  retaining deterministic omics-file order and subsetting both inputs to the
  same analysis sample set.
- The console and run log report aligned, omics-only, and phenotype-only sample
  counts. Inputs with no sample IDs in common still fail explicitly.
- Added end-to-end CLI coverage for unequal ID sets and verified the behavior
  against a 5,100-sample TOPMed expression matrix.

## Statistical and modeling extensions

- Added complete-case mixed-formula compilation, sparse pedigree combinations,
  richer covariance profiles, joint Kenward-Roger tests, and expanded
  multilevel and sparse meta-analysis.
- Added bounded multidimensional adaptive Gauss-Hermite quadrature and
  deterministic multi-start diagnostics for sparse ZIP/ZINB fits.
- Added covariance-aware multivariable MR conditional-strength inference and
  strengthened positive-mixture and SKAT-O tail calibration.

## Omics, runtime, and validation

- Includes the v0.3.1 R-compatible `winsor_mad(k=4)` transform and its complete
  CLI tutorial, with whitespace-tolerant transform pipeline syntax.
- Upgraded the pinned runtime to JDistlib 0.10.2 with checksum verification and
  independent downstream numerical validation.
- Expanded Java/R fixtures, boundary tests, GPU runtime guidance, and
  reproducible audit documentation.

Release downloads: `jlinalg-0.3.2.jar` is the self-contained executable;
`JLinAlg-0.3.2-library.jar` is the thin library. Sources and Javadoc JARs and
`SHA256SUMS.txt` are also included.

# JLinAlg 0.3.0

## Post-release development

- Mixed formulas can map `(1|individual)` grouping columns to pedigree
  precision terms, retain unobserved pedigree members, and combine them with
  ordinary independent random effects. Optional complete-case compilation
  applies one row mask to fixed/random predictors, grouping columns, weights,
  and offsets and exposes the retained original rows.
- Correlation likelihood profiles now accept any selected effect pair in
  unstructured random blocks larger than two while retaining exact plus/minus
  one boundary coordinates.
- Sparse meta-analysis accepts a known sparse sampling precision and estimates
  independent REML or Paule-Mandel heterogeneity without forming its inverse.
  Multilevel meta-analysis adds compound-symmetry and AR(1) covariance families
  and nuisance-refitted random-standard-deviation profiles.
- Sparse mixed models add joint Kenward-Roger F tests, one-boundary mixture
  profile intervals, original-coordinate within-group PEV blocks, and joint
  unstructured grouped plus sparse-pedigree covariance optimization.
- Multivariable MR adds covariance-aware generalized estimation and true
  conditional-strength statistics when per-instrument cross-exposure
  covariance is supplied.
- Added bounded multidimensional adaptive Gauss-Hermite quadrature for
  binomial-logit and Poisson-log crossed or sparse-precision random effects.
- Sparse ZIP/ZINB fits now use deterministic multi-start outer optimization and
  report start count, distinct stationary modes, and whether the selected mode
  was reproduced and best among the deterministic starts. Fixed-effect
  covariance conditions on nuisance variances at active bounds.
- Positive-mixture tails fall back to a labeled saddlepoint calculation when
  the bounded gamma series is unavailable. SKAT-O now defaults to parametric
  Gaussian-score-null simulation; analytic moment matching is opt-in.

## Remaining-code accuracy and performance audit

- Corrected shared statistical scaling, exact Kendall tails, rank-test cutoffs,
  zero-sign ties, unavailable inference, and LOESS weight/leverage normalization.
- Corrected Gaussian GLM likelihood and finite-DF inference, unit-sensitive rank
  and penalty handling, GEE variance scaling, and cross-validation accounting.
- Replaced conditional-beta/Fisher-curvature Laplace fitting with full marginal
  optimization and observed-mode curvature; corrected beta mixed-model likelihood
  agreement and rejected unidentifiable REML covariance directions.
- Repaired Cox risk-set cancellation, nonlinear convergence/final linearization,
  sparse assembly, and scale-safe mediation delta inference.
- Repaired colocalization support and weighted evidence, SuSiE variance state,
  tiny conjugate means, LD/GRM validity checks, MR centering/scaling, and unphased
  LD global stationary-point selection without multi-start timeouts.
- Repaired streaming failure/row accounting, annotation keys, dosage validation,
  bounded multi-pass BH merging, and quadratic-form tail scaling.
- Added independent R/reference fixtures, paired warm timing reports, new website
  guides for basic tests, nonlinear models, mediation, and beta regression, and
  versioned library/source/Javadoc/CLI release assets with SHA-256 checksums.

Compatibility correction: `MultivariableMrResult.conditionalFStatistics()`
previously returned marginal mean z², not conditional strength. It is deprecated
and now fails explicitly. Use `marginalFStatistics()` for the actual statistic;
covariance-aware fits now populate it with true conditional strength, also
available with its Q statistics from `conditionalStrength()`.

The [release audit](https://github.com/robbyjo/JLinAlg/blob/v0.3.0/docs/release-0.3.0-audit.md)
records reproduction commands, measured faster/slower cases, test gates, and
retained approximation limits. Earlier timing claims from incorrectly optimized
likelihoods are not performance baselines for the repaired methods.

Release downloads: `jlinalg-0.3.0.jar` is the self-contained executable;
`JLinAlg-0.3.0-library.jar` is the thin library. Sources and Javadoc JARs and
`SHA256SUMS.txt` are also included. The library classifier prevents a
case-insensitive filename collision on Windows.

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
  [the validation report](https://github.com/robbyjo/JLinAlg/blob/v0.3.0/docs/advanced-validation.md).

## Earlier additions

- Added sparse first-order Laplace beta mixed models with analytic beta
  scores and observed curvature, jointly optimized fixed effects, precision and variance
  components, grouped random effects, arbitrary sparse coefficient precision,
  and a pedigree convenience API consuming sparse `A^-1`. An independent
  `glmmTMB` 1.1.14 fixture checks fixed effects, precision, variance, likelihood,
  unrelated-founder equivalence, and retention of unphenotyped ancestors.
  Corrected marginal optimization reduces the original grouped log-likelihood
  error from about 0.122 to below 1e-8; the release audit supersedes prior
  timings from the incorrect objective.
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
