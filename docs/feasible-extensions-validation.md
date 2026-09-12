# Quantile, missing ordinal, and diffuse ARIMA extensions

Validated 2026-09-12. This first batch implements three bounded extensions from
the open-method inventory. Larger estimator changes remain in [TODO](../TODO.md);
no global-optimality, bandwidth-validity, or finite-sample calibration limit is
removed by these changes.

## Implemented contracts

| API | Added behavior | Remaining boundary |
| --- | --- | --- |
| `QuantileRegressionInference` | Exact LP fit with full asymptotic covariance from supplied conditional densities; optional pooled iid Gaussian residual kernel at caller-set bandwidth | Independent, continuous responses and regular positive densities; no automatic inference-valid bandwidth, dependent-row or mass-point inference |
| `SemOrdinal.fitPairwiseMissing` | `-1` missing categories, available-pair PML, joint threshold/structural case or cluster covariance | MCAR or justified pair-specific observation mechanism; no general MAR ordinal FIML, mixed continuous/ordinal estimator, DWLS/WLSMV, multigroup or ordinal modification indices |
| `DiffuseArima.Result` | Full observed-information coefficient covariance and SEs, including seasonal terms and drift/dynamic cross terms | Gaussian-likelihood asymptotics; no general regressors or historical smoothing; unresolved information and transform-bound fits suppress inference |

Quantile inference uses a separate result wrapper so the exact LP result and its
optimality certificate retain their existing meaning. Nonconverged exact fits
cannot enter the new covariance API. Both ordinal fitting paths now suppress
all covariance and parameter inference on nonconvergence.

Missing ordinal rows with fewer than two observed responses contribute no pairs
and are removed. Cluster IDs index the original input, then follow retained
rows. Every pair must have at least two jointly observed rows. A declared
category must occur in the retained sample. Complete data remain accepted by the
new API; the existing complete-data API continues to reject missing categories.

Diffuse inference uses the same exact diffuse likelihood as fitting. Drift is
profiled for optimization but included jointly for covariance, so its uncertainty
and cross terms are retained. Innovation variance remains profiled out. Drift
location and its SE have differenced-process units; divide by the seasonal
period for slope-per-time-step units when seasonal differencing is used. Two
finite-difference step sizes test information stability, with no ridge repair.
The old four-argument result constructor is retained and marks covariance
unavailable. The original optimizer evaluation count excludes Hessian work.

## Independent numerical evidence

Fixtures were generated locally using R 4.6.1, lavaan 0.7.2, and quantreg 6.1.
The generator and small frozen data/results are checked in; Java tests do not
require R or package downloads.

- **Quantile:** independent `quantreg::rq(method="fn")` coefficients and R
  matrix-algebra density sandwiches, with varying densities and a pooled kernel.
  Tests also check the analytic normal-median variance `pi/(2*n)`, joint response
  and predictor unit transformations, defensive result copies, and rejection
  of invalid densities, deficient designs, and nonconverged exact fits.
- **Ordinal:** lavaan PML with `missing="pairwise"` on 588 retained rows of a
  frozen four-indicator ordinal sample. Independent balanced binary tables with
  disjoint observed pairs check the analytic probit correlation and its variance.
  Tests verify ignored empty/single-response rows, cluster replication invariance,
  missing-code/coverage validation, and the complete-data API's strict behavior.
- **Diffuse:** nine `stats::arima` covariance fixtures cover AR, MA, ARMA,
  seasonal AR/MA, twice-integrated, missing, ordinary drift, seasonal drift, and
  monthly missing series. Every covariance entry must agree within 0.002 times
  the geometric mean of its reference diagonal variances. This allows R's finite
  diffuse initialization and independent numerical Hessian differences.
  Analytic irregular-gap random-walk GLS drift variance, response scaling,
  compatibility construction, empty coefficient sets, and nonconvergence are
  checked separately.

The quantile covariance formula follows the sandwich described in the
[quantreg reference manual](https://stat.ethz.ch/CRAN/web/packages/quantreg/quantreg.pdf).
The ordinal estimator is an available-pair composite likelihood; lavaan provides
[pairwise-likelihood tutorials](https://lavaan.ugent.be/resources/tutorial.html).
Its case/cluster sandwich does not supply WLSMV statistics or an ordinal full
likelihood. See the [SEM](vignettes/sem.md), [time-series](vignettes/time-series.md),
and [quantile](vignettes/regression-families.md) vignettes for API examples.

## Measured full-fit workloads

Local Windows run on an Intel Core i9-14900K, three warmups and seven measurements. Values below are
medians for complete fitting plus covariance. Input parsing and accuracy gates
are outside the measured region. These are workload measurements, not claims of
speed relative to R or previous implementations.

| Workload | Median ms | Maximum absolute covariance error vs R |
| --- | ---: | ---: |
| Quantile, supplied densities, 101 rows | 0.790200 | 8.33e-17 |
| Quantile, iid kernel, 101 rows | 0.486500 | 2.64e-12 |
| Missing ordinal PML, 588 rows, 4 indicators | 5.518800 | 4.43e-8 |
| Missing diffuse AR(1) with drift, 240 positions | 2.423400 | 3.29e-6 |

Quantile covariance stores O(np + p²) arrays. Diffuse information adds O(k²)
filter evaluations while keeping each covariance workspace state-sized.
Ordinal inference retains the existing dense parameter Hessian and per-case
score storage; this batch makes no new high-dimensional performance claim.

## Reproduction

From the repository root:

```powershell
# Optional regeneration; requires the named R packages already installed.
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/test/resources/r-reference/generate-feasible-extensions.R

.\gradlew.bat check benchmarkClasses
.\gradlew.bat benchmarkFeasibleExtensions
.\gradlew.bat javadoc
```

The new tests are `QuantileRegressionInferenceTest`, `SemMissingOrdinalTest`,
and `DiffuseArimaInferenceTest`. The final full gate ran 719 tests: 716 passed,
zero failures/errors, and three optional native CHOLMOD skips. Website checks,
benchmark compilation, Javadoc generation, and the new accuracy-gated benchmark
also passed.
