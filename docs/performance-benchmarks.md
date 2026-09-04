# Performance benchmarks

`benchmarkAssociation` runs deterministic end-to-end macrobenchmarks for:

- covariate-reused block OLS;
- prepared-null block GLM score association;
- exact parallel per-predictor OLS refits.

`benchmarkMixedModels` measures sparse-equation REML and reports both equation
and factor nonzero counts. For workloads with at most 2,000 observations it
also runs the dense covariance reference path.

The task performs warm-up iterations, reports the median of measured runs, and
consumes a result coefficient to prevent dead-code elimination. It is intended
for workload and backend comparisons rather than nanosecond-scale JVM method
measurement.

```powershell
.\gradlew.bat benchmarkAssociation
.\gradlew.bat benchmarkMixedModels
.\gradlew.bat benchmarkTimeSeries
.\gradlew.bat benchmarkXwasMr
.\gradlew.bat benchmarkTopmedPenalized
```

The real-data ridge/LASSO/elastic-net comparison with R `glmnet` is documented
in [TOPMed penalized-regression performance](topmed-penalized-performance.md).
The bounded parallel MR screen and its independent R numerical reference are
documented in [xWAS MR CLI and performance](xwas-mr-cli-performance.md).

## Time-series benchmarks

`benchmarkTimeSeries` loads six fixed real-world series from
`src/benchmark/resources/timeseries`. The cases cover conditional AR, ARMA,
ARIMA, seasonal ARIMA, exhaustive small-order selection, and exact stationary
AR likelihood. Resource headers, row order, finite values, row counts, and
SHA-256 checksums are validated before timing begins.

Run a comma-separated subset while profiling a particular implementation:

```powershell
.\gradlew.bat `
  "-Djlinalg.benchmark.cases=conditional_arma21_sunspots,exact_ar2_nile" `
  "-Djlinalg.benchmark.warmups=3" `
  "-Djlinalg.benchmark.measurements=7" `
  benchmarkTimeSeries
```

Available cases are:

- `conditional_ar2_nile`
- `conditional_arima310_wwwusage`
- `conditional_sarima_airline`
- `conditional_sarima_ukgas`
- `conditional_seasonal_ar_nottem`
- `conditional_arma21_sunspots`
- `automatic_nile`
- `exact_ar2_nile`

Output is CSV with median wall time, optimizer evaluations, and convergence.
The automatic-order case reports `-1` evaluations because its public candidate
results do not expose a total. Dataset provenance, source links, checksums, and
licensing notes are in the
[time-series data README](../src/benchmark/resources/timeseries/README.md).

For an R `stats::arima` comparison, run the bundled script. Conditional cases
use `method = "CSS"`; `exact_ar2_nile` uses `method = "ML"`. Its arguments are
data directory, warm-ups, measurement batches, repetitions per batch, and an
optional comma-separated case filter:

```powershell
Rscript src/benchmark/r/time_series_benchmark.R `
  src/benchmark/resources/timeseries 1 7 50
```

Repeated fits are used because Windows R elapsed-time resolution is too coarse
for the shortest cases. Compare `median_seconds` from Java with
`median_seconds_per_fit` from R on the same otherwise-idle machine.

Dimensions and execution can be controlled without editing source:

```powershell
.\gradlew.bat `
  "-Djlinalg.benchmark.rows=50000" `
  "-Djlinalg.benchmark.variables=10000" `
  "-Djlinalg.benchmark.covariates=10" `
  "-Djlinalg.benchmark.parallelism=16" `
  "-Djlinalg.benchmark.warmups=3" `
  "-Djlinalg.benchmark.measurements=7" `
  benchmarkAssociation
```

Output is CSV with median wall time and variables per second. Record the Java
version, CPU/GPU model, backend, thread environment, heap settings, and power
mode alongside published results. Run on an otherwise idle machine and compare
multiple repetitions before changing routing thresholds.

The default workload deliberately fits only 100 exact OLS models because that
path is a correctness/per-fit baseline; the fast paths scan the full requested
variable count.
