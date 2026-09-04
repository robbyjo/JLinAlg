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
.\gradlew.bat benchmarkColocSusie
.\gradlew.bat benchmarkSusie
.\gradlew.bat benchmarkTopmedSem
```

The real-data ridge/LASSO/elastic-net comparison with R `glmnet` is documented
in [TOPMed penalized-regression performance](topmed-penalized-performance.md).
The bounded parallel MR screen and its independent R numerical reference are
documented in [xWAS MR CLI and performance](xwas-mr-cli-performance.md).
The observed-variable path-model comparison with R `lavaan` is documented in
[TOPMed SEM validation and performance](topmed-sem-performance.md).

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

## SuSiE benchmark

`benchmarkSusie` runs the complete individual-data fit on response replicate 1
of the 574-by-1,001 `N3finemapping` data from the official susieR fine-mapping
vignette. The timed region includes centering, scaling, sufficient-statistic
construction, nine IBSS iterations, and credible-set calculation; data loading
is excluded. The matching R script uses the same fixed absolute prior variance,
residual-variance estimation, tolerance, and credible-set settings:

```powershell
.\gradlew.bat benchmarkSusie
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/benchmark/r/susie_vignette_benchmark.R
```

On Windows 11 with an Intel Core i9-9900K, Eclipse Adoptium Java 21.0.4, and
R 4.6.1/susieR 0.14.2, three warmups followed by seven measurements produced
median end-to-end times of 0.0704391 seconds for JLinAlg's portable CPU path
and 0.340 seconds for R: JLinAlg was 4.83 times faster. Both fits converged in
nine iterations, returned three credible sets, and estimated residual variance
as 6.36378419836.

The speedup comes from computing only one triangle of the symmetric cross-product
in parallel, caching every effect's `X'X b` contribution, and reusing it for
both residual and ELBO updates. Override warmups, measurements, or backend with
`jlinalg.benchmark.susie.warmups`, `jlinalg.benchmark.susie.measurements`, and
`jlinalg.benchmark.susie.backend`.

## SuSiE colocalization benchmark

`benchmarkColocSusie` measures the Bayes-factor combination step after SuSiE
fine mapping. Its default workload uses 10,000 variants and 10 signals per
trait (100 signal pairs), with deterministic inputs, three warmups, and the
median of seven measured runs. Override dimensions with
`jlinalg.benchmark.coloc.variants` and `jlinalg.benchmark.coloc.signals`.

On an Intel Core i9-9900K with Java 25, two verification runs completed in
median times of 0.0236 and 0.0279 seconds, or 35.9-42.3 million
signal-pair/variant updates per second.
The implementation performs variant matching and per-signal log-sums once,
then uses contiguous primitive arrays in the pairwise hot loop. The returned
100-by-10,000 conditional-H4 posterior matrix is included in the timing.

The default workload deliberately fits only 100 exact OLS models because that
path is a correctness/per-fit baseline; the fast paths scan the full requested
variable count.
