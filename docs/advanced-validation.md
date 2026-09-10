# Statistical audit: accuracy, performance, and remaining limits

Run date: 2026-09-08. This report supersedes the September 7 descriptions that
treated partial implementations as completed statistical methods. The repair
set includes production algorithms, independent R fixtures, regression tests,
consumed-output benchmarks, source vignettes, and the static website.

## What was repaired

| Area | Delivered method and independent evidence |
|---|---|
| SEM | Joint latent RAM and means; constrained observed-pattern FIML; ordinal probit PML; full covariance, robust/cluster inference, efficient continuous-model modification indices, indirect delta inference. Joint lavaan estimate error below 7.4e-8, covariance below 3.8e-9, likelihood below 9e-12 on the fixtures. |
| Meta-analysis | Joint ML/REML random-moderator covariance, correlated CR0/CR1/CR2, finite DF, PET/PEESE, iterative L0/R0/Q0 trim-and-fill, stable Hedges correction. 236 metafor/clubSandwich reference values; random covariance error below 7e-7 and likelihood below 1e-8. |
| Mixed models | Multiple sparse unstructured blocks, model-derived Satterthwaite/KR derivatives, weights/offsets, constrained formula-native ML profiles. Four R fixtures: likelihood below 2e-11, covariance entries below 0.00085, SE below 4e-6, finite DF below 0.0006; profile endpoints within 0.0001. |
| Time series | Symbolic exact diffuse Kalman filtering for integrated/seasonal/missing models; state-sized missing-series likelihood and forecasts. Thirteen stats::arima fixtures: coefficients within 2.38e-6, forecast means within 6.23e-6, forecast SE within 1.27e-5. |
| MR | Actual conditional instrument association and forward/backward signal selection, generalized/overlap-aware CLI, confidence propagation, converged-iterate uncertainty, and corrected SVG Egger line. Independent R conditional/GLS calculations agree within 2e-11. |
| GLMM | Mode-adaptive, node-refined scalar random-intercept integration, trials/offsets, zero-variance comparison, joint observed covariance. Eight R integrate likelihoods within 1e-12; four lme4 AGQ fits have coefficient/SD errors below 4e-7 and SE errors below 2e-7. |
| Regression | Shared multivariate QR, corrected smoothed quantile gradient plus certified nonsmoothed LP, score-checked multinomial optimizer, stable kernel weights, R-compatible supersmoother, smoother-aware Robinson slope inference. 1,006-value R fixture also covers selection inference; 25 additional exact-quantile cases have objective error below 1.35e-9. |
| Penalized inference | Weighted independent held-out inference plus fixed-penalty Gaussian active-set/sign polyhedral inference for LASSO/elastic net. glmnet/selectiveInference LASSO references, analytic elastic-net truncation, null calibration, and extreme-tail tests. |

These are errors on named fixtures, not global accuracy bounds. Optimizer
convergence is required but does not establish identification, global optimality,
or substantive model fit. The source vignettes specify supported assumptions.

The final integrated gate discovered **462 tests: 459 passed, zero failed,
three optional native CHOLMOD tests skipped**. All benchmark classes and the
executable CLI JAR compiled; website structure/link/theme checks passed. The
earlier 406-test green suite was insufficient: independent review found the
counterexamples below, which were repaired and added as regression checks.
Benchmark timing claims use separately recorded runs, not concurrent Gradle
tasks. Native CHOLMOD execution is not certified by this portable run.

### Counterexamples added after independent review

- SEM: near-unit ordinal correlations and probabilities down to 7.36e-74,
  underidentified paths/latent means, empty ordinal cells, mean-inclusive SRMR,
  large-location FIML moments, invalid auxiliary saturated diagnostics, and
  strict backend policies.
- Mixed: adding 1e8 to a response preserves variances/likelihood; configured
  physical variance bounds are enforced; correlations at plus/minus one have
  available-side profiles. Categorical slopes have 40 design comparisons and
  five weighted/offset lme4 fit comparisons.
- Time series: leading missing prefixes up to 100,000 positions, unrestricted
  drift with sparse observations, multimodal scalar likelihoods, very small
  response units, optional initializer failure, unidentified AR/MA cancellation,
  and rejection of infinite observations.
- GLMM: concentrated-binomial nonzero variance near the boundary, stable
  count probabilities through 1e16, ULP-aware count validation, and Poisson
  posterior modes with offsets up to 1000.
- Regression/inference: scaled selection constraints, near-boundary selective
  intervals, tiny representable p-values, smoothing-induced noise covariance,
  kernel weighted-sum overflow, and the periodic smoother minimum window.

The exact quantile API adds primal/dual feasibility, complementarity, and
duality-gap checks for the nonsmoothed LP; nonunique solutions are assessed by
objective and optimality rather than demanding one arbitrary coefficient vector.

## Measured speed, including slower cases

Windows development host, Java 25 portable CPU, R 4.6.1. Inputs are identical
within each comparison and prepared outside the measured fit. Warmup policy,
replicates, checksums, package versions, and estimator caveats are recorded in
the linked family reports. The host was shared, not a controlled performance
laboratory. Public APIs may perform different ancillary bookkeeping.

| Workload | Java | R | Interpretation |
|---|---:|---:|---|
| Multivariate OLS, 512 rows | 0.106 ms | 0.244 ms | 2.30x R/Java on this fit |
| Multinomial, 512 rows | 0.605 ms | 2.867 ms | 4.74x R/Java |
| Matching smoothed quantile, 512 rows | 0.294 ms | 0.993 ms | 3.38x versus independent smoothed R optim, not nonsmoothed rq |
| Exact nonsmoothed quantile, 1024 rows / 6 columns | 1.818 ms | 0.631 ms | R Frisch–Newton faster; Java includes original-unit optimality certificate |
| Fixed-bandwidth kernel, 512 rows | 1.443 ms | 7.753 ms | npreg includes its usual result construction |
| Supersmoother, 512 rows | 0.106 ms | 0.098 ms | R approximately 7% faster |
| Robinson + HC3, 512 rows | 8.200 ms | 13.047 ms | Compared to explicit R matrix implementation |
| Latent SEM fixture | 0.957 ms | 26 ms | lavaan includes higher-level model bookkeeping |
| FIML SEM fixture | 6.494 ms | 112 ms | Same constrained likelihood |
| Random-moderator REML, 36 effects | 4.44 ms | 37.3 ms | Same ML/REML estimand; small dense workload |
| Correlated GLS + CR2 | 0.210 ms | 8.50 ms | Coefficient covariance and t inference |
| Sparse mixed Satterthwaite, 12,000 rows | approximately 0.85 s | 0.220 s | Java slower; Java fits under a 96 MB heap |
| Bernoulli random-intercept AGQ | 111.397 ms | 80 ms | Java slower; stricter refinement/multistart controls |
| Rare Bernoulli AGQ | 208.816 ms | 110 ms | Java slower |
| Poisson random-intercept AGQ | 69.338 ms | 110 ms | Java faster on this fixture |
| Stationary ARMA likelihood, 240 positions | 0.0108 ms | 0.0253 ms | Fixed coefficients, filtering only |
| Missing monthly seasonal likelihood | 0.3974 ms | 0.3429 ms | Java about 16% slower |
| Conditional associations, 20,000 five-SNP fits | 0.060334 s total | 2.08 s total | Preconstructed summary-score calculations |
| Fixed-penalty selective inference, 200 n=160 fits | 0.024986 s total | 0.510 s total | Tight CI inversion; R first computes its native coarse grid |

Java is not uniformly faster. Sparse factorization reuse reduced the large
mixed-model implementation from roughly 13 seconds to under one second,
but R's Satterthwaite fit is still faster. The categorical mixed-model fits
measured 0.07–2.64 seconds in Java versus 0.01–0.19 seconds in R. No claim is made for speed outside
these workloads or for end-to-end file/CLI startup throughput.

Two important comparison differences must not be hidden:

- R ARIMA uses finite diffuse initialization, while Java uses symbolic exact
  diffuse recursion. The maximum likelihood difference is 8.67e-5 at R
  kappa=1e6 and drops to 8.68e-6 at kappa=1e7. Diffuse full-fit R timings include
  coefficient covariance that the Java diffuse API does not expose; use the
  fixed-parameter likelihood timings for the narrower kernel comparison.
- Java's quantile objective is smoothed. Nonsmoothed rq took 0.387 ms on the
  regression fixture but is not the denominator in the smoothed-loss ratio.
  R selectiveInference uses a coarse CI search; the fixture additionally
  inverts its reported truncation distribution with R uniroot for accurate
  endpoint comparison.

## Evidence and reproduction

- [Regression and selection accuracy/timing](../src/benchmark/resources/regression/accuracy-and-timing.md)
- [Exact nonsmoothed quantile LP accuracy/timing](../src/benchmark/resources/exact-quantile-benchmark/exact-quantile-results.md)
- [Selective inference boundaries and paired timing](../src/benchmark/resources/r-reference/selective-inference-notes.md)
- [Joint SEM references, timing, and assumptions](vignettes/sem.md)
- [Meta-analysis references and paired measurements](vignettes/meta-analysis.md)
- [Sparse mixed-model algorithms, R profiles, timings, memory limits](lme4-pedigreemm-roadmap.md)
- [Exact ARIMA measurements and reproduction](../src/benchmark/resources/timeseries/exact-arima-measurements.md)
- [Conditional MR assumptions, fixtures, and timing](vignettes/mr-end-to-end.md)
- [Adaptive quadrature validation and timing](../src/benchmark/resources/quadrature/accuracy-and-timing.md)
- [Feature contracts and API navigation](vignettes/advanced-extensions.md)

```powershell
./gradlew.bat check benchmarkClasses executableJar
# Run one benchmark at a time for measurements:
./gradlew.bat benchmarkRegressionFamilies
./gradlew.bat benchmarkExactQuantile
./gradlew.bat benchmarkSemJoint
./gradlew.bat benchmarkMetaParity
./gradlew.bat benchmarkSparseMixedParity
./gradlew.bat benchmarkMrConditional
./gradlew.bat benchmarkSelectiveInference
./gradlew.bat benchmarkGlmmQuadrature
./gradlew.bat benchmarkExactTimeSeries
```

The family reports link the exact R generators and runners. Tests use frozen
fixtures and do not require R or network access. Regeneration requires the
listed R packages; scripts use the optional repository-local build/r-library
without replacing the user's library. Neither startup nor package installation
is included in the benchmark figures.

## Remaining scope

See [TODO.md](../TODO.md) for the explicit remaining extensions. Important
limits include ordinal DWLS/missing-ordinal/multigroup SEM; sparse meta models
beyond a known precision plus independent heterogeneity; higher-order mixed
variance-boundary intersections; diffuse coefficient covariance; non-Gaussian
conditional-GWAS summary models; tensor growth in multidimensional AGQ;
single-mode conditional Laplace integration for zero-inflated fits;
unknown-noise or response-tuned polyhedral inference; and quantile inferential
covariance/mixed-type kernel models.
The completed implementations must not be advertised as full package parity.
