# Exact ARIMA validation and timing — 2026-09-08

Measured on Windows, Intel Core i9-14900K, Oracle JDK 25, R 4.6.1, CPU filtering.
The shared host was not otherwise isolated; these are local wall-clock
measurements. All 13 series have 240 original positions. Inputs and R-generated
expected values are frozen under `src/test/resources/timeseries`.

The accompanying `exact-arima-measurements.csv` contains all measured timings
and fit errors. No conditional-likelihood timing is mixed into these results.
These measurements were refreshed after repairing all seven independent-review
findings, including unrestricted GLS drift profiling and scalar basin search.

| Fixture | Java likelihood (ms) | R likelihood (ms) | Java full fit (ms) | R full fit (ms) |
|---|---:|---:|---:|---:|
| Stationary ARMA(1,1) | 0.0108 | 0.0253 | 4.177 | 1.552 |
| Missing stationary ARMA(1,1) | 0.0110 | 0.0247 | 4.519 | 1.521 |
| ARIMA(1,1,0) | 0.0117 | 0.0172 | 2.106 | 0.936 |
| ARIMA(0,1,1) | 0.0337 | 0.0264 | 1.807 | 0.782 |
| Missing ARIMA(1,1,1) | 0.0188 | 0.0518 | 1.243 | 2.058 |
| Seasonal (0,1,1)(0,1,1)[4] | 0.0866 | 0.0829 | 4.615 | 6.567 |
| Missing seasonal | 0.0881 | 0.0832 | 3.702 | 5.396 |
| Seasonal-only (0,0,0)(1,1,0)[4] | 0.0594 | 0.0591 | 5.729 | 2.216 |
| Missing ARIMA(0,2,1) | 0.0222 | 0.0301 | 2.219 | 1.355 |
| Missing ordinary drift | 0.0120 | 0.0183 | 2.113 | 2.676 |
| Missing seasonal drift | 0.0545 | 0.0617 | 6.546 | 7.123 |
| Missing (1,1,1)(1,1,1)[4] | 0.1094 | 0.0937 | 16.700 | 23.090 |
| Missing monthly (0,1,1)(0,1,1)[12] | 0.3974 | 0.3429 | 18.115 | 26.198 |

Fixed-parameter likelihood measurements include state initialization and
filtering, with polynomial expansion and input loading outside the clock. They
exclude optimization, coefficient Hessians and forecasts. Each median is over
11 batches of 1,000 evaluations following three warmup batches; every likelihood
is consumed by a checksum. R uses the installed `stats:::makeARIMA` and
`C_ARIMA_Like`, with Rossignol stationary initialization and `kappa=1e6`.
Java uses matrix doubling for the stationary block and symbolic exact diffuse
recursions. Both evaluate the same Gaussian estimand at the frozen R
coefficients, subject to R's finite-diffuse approximation. R's drift series is
centered before the clock; Java subtracts the fixed trend during filtering.

These results do not establish universal superiority over R. The Java kernel is
about 1.1–2.8 times faster on several small/nonseasonal cases, approximately tied
on some small seasonal models, and about 16% slower on the monthly missing case.
The integrated MA kernel was about 28% slower in this run. Shared-host timing
variation is visible, so compare the full table rather than a single speed ratio.
Dense covariance updates are still state-sized; their cost grows with model
order and seasonal period.

Full-fit medians use ten warmups and 31 measured runs. All fitted parameter
checksums are consumed; all fits reported convergence and passed accuracy
gates. Integrated fits include eight-step forecasts. Java's stationary fits
use five starts, integrated fits one start, and tolerances 1e-8 and 1e-9,
respectively. R uses BFGS, one start, `reltol=1e-12`, and `maxit=2000`.
One-dimensional Java objectives now explore 64 intervals and refine sampled
basins instead of assuming unimodality. This increases full-fit time for several
one-parameter models. Drift is now profiled analytically, reducing the nonlinear
dimension while preserving its likelihood parameter count. Stationary covariance
validation evaluates information at two step sizes without ridge regularization.
**Full-fit inference work is different**: R calculates coefficient covariance;
Java does so for stationary fits but does not yet expose this for diffuse fits.
These full-fit timings therefore describe the APIs as delivered and are not a
claim of equal optimizer or inference workload.

Maximum optimized coefficient error was 2.38e-6; maximum forecast mean and SE
errors were 6.23e-6 and 1.27e-5. The largest ML log-likelihood difference was
8.67e-5. At identical frozen coefficients and R `kappa=1e7`, the maximum
likelihood difference decreases to 8.68e-6. This supports the finite-diffuse
approximation explanation; Java's exact-diffuse normalization and missing-data
semantics are explained in `docs/vignettes/time-series.md`.

All 35 time-series JUnit tests passed with a 128 MB Java heap. Tests include
fixed-coefficient R likelihood and forecast comparisons, optimized R fits,
an analytic random-walk missing-gap likelihood, a dense stationary oracle,
complete-difference equivalence, near-unit stationary initialization, drift
without any observable adjacent differences, and a 50,000-position missing
series. The latter would require 20 GB for the old original-series covariance
alone and now fits within the 128 MB test process. Eight additional tests cover
the seven independent-review findings: invariance to 100,000 leading NaNs,
ordinary and seasonal sparse-observation drift MLEs of 100 and 400, the competing
MA(1) optimum near -0.160012, positive variances at tiny scales, optional
conditional initialization for short exact AR(2), unavailable inference for
unidentified AR/MA cancellation, and rejection of either sign of infinity.
The additional R generator records exact likelihood references for the small
MA/AR examples and the analytic drift MLEs. No project-wide build was
run while parallel packages were being edited.

## Reproduction (PowerShell, repository root)

The R scripts begin with `.libPaths(c('build/r-library', .libPaths()))` and need
only base `stats`. Full fixture generation is deterministic for this R version.

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/timeseries/generate-exact-arima.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/timeseries/generate-review-arima.R

$tsCp = 'build/timeseries-classes;build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar'
& 'C:/Program Files/Java/jdk-25/bin/javac.exe' -cp $tsCp -d build/timeseries-classes (Get-ChildItem src/main/java/org/jlinalg/timeseries/*.java).FullName src/benchmark/java/org/jlinalg/benchmark/ExactTimeSeriesBenchmark.java src/benchmark/java/org/jlinalg/timeseries/ExactArimaLikelihoodBenchmark.java

java -cp $tsCp org.jlinalg.timeseries.ExactArimaLikelihoodBenchmark
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/exact_arima_likelihood_benchmark.R
java '-Djlinalg.benchmark.warmups=10' '-Djlinalg.benchmark.measurements=31' -cp $tsCp org.jlinalg.benchmark.ExactTimeSeriesBenchmark
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/exact_time_series_benchmark.R 10 31

# Once all parallel package edits are stable:
./gradlew.bat test --tests 'org.jlinalg.timeseries.*'
```

Isolated test execution used local copies of the cached JUnit 5.13.4 and
JUnit Platform 1.13.4 jars, `javac` output under `build/timeseries-test-classes`,
and a launcher selecting `org.jlinalg.timeseries`. Test resources were supplied
from `src/test/resources`. Native backend tests emitted existing JDistlib GPU
initialization notices; the likelihood and benchmark loops themselves run on
CPU. Final full-project Gradle verification is left to the coordinating parent.
