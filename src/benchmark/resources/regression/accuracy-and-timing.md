# Regression repair validation — 2026-09-08

Windows, portable Java CPU backend, Oracle JDK 25, R 4.6.1. Reference packages:
nnet 7.3-20, quantreg 6.1, np 0.70-5, glmnet 5.0, selectiveInference 1.2.5.
The shared development host was not isolated from other work. These are local
wall-clock observations, not universal performance promises.

## Accuracy gates

`RegressionAccuracyTest`, `RegressionRParityTest`, `PartiallyLinearNoiseTest`, and `SelectionAccuracyTest`
passed all 19 isolated tests before full-project integration. Their frozen R
fixture contains 1,006 values. The comparisons include:

- Multivariate coefficients and residual covariance: absolute tolerance 1e-12.
- Multinomial coefficients: 1e-5; full log likelihood: 1e-9; score convergence required.
- Quantile coefficients: 1e-5 against independently optimized identical smoothed
  loss; 0.003 against nonsmoothed `quantreg::rq` for smoothing 0.001. These are
  distinct objectives. The original median `[1,2,3]` failure now returns 2.
- Gaussian kernel fitted values: 1e-12 against direct R kernel algebra.
- Supersmoother: 1e-10 against `stats::supsmu`, covering default, weighted/bass,
  fixed span and periodic curves. The original 40-row example matches to the
  printed double precision; repeated queries and tied extrapolation are tested.
- Robinson partialling-out coefficients, centered smooth, smoother-aware HC3-style
  and homoskedastic covariance: 1e-12 against independent R matrix algebra.
  Original-observation influence and smoothing-induced noise correlations are
  included. An analytic paired-noise test checks the homoskedastic sampling
  variance; generalized HC3 is conservative at high smoothing leverage.
- Conditional LASSO: R `fixedLassoInf` coefficient targets, truncation limits and
  conditional confidence inversion; two-sided p-values are compared after
  converting R's sign-oriented one-sided p-values. R's default coarse confidence
  grid misses its target tails by up to about 0.0009 on this fixture, so the R
  generator additionally inverts its reported truncation distribution with
  `uniroot(tol=1e-12)`. Java matches the refined intervals within 1e-5 and
  truncation limits within 1e-6.
- Analytic LASSO/elastic-net truncation, conditional p-value calibration,
  extreme-tail nonzero p-values, weighted held-out inference, and p>n selection.

These are explicit test tolerances, not claims of worst-case errors over all
possible datasets. Identifiability, bandwidth assumptions, fixed-lambda/known
Gaussian noise requirements and optimizer convergence still matter.

## Measured fit timings

Identical deterministic inputs, n=512, two design columns and two OLS outcomes.
Java: 100 warmup fits and median of 101 measured fits per operation. R: five
warmups, seven batches of 100 fits, median per-fit wall time. Startup and input
construction are excluded. Each timed result is consumed in a checksum;
iterative fits must pass convergence checks. Java and R were run sequentially.

| Estimator | Java ms | R ms | R / Java |
|---|---:|---:|---:|
| Multivariate OLS / `lm` | 0.106 | 0.244 | 2.30 |
| Multinomial / `nnet::multinom` | 0.605 | 2.867 | 4.74 |
| Smoothed quantile / R `optim` | 0.294 | 0.993 | 3.38 |
| Gaussian kernel / `np::npreg` with fixed bandwidth | 1.443 | 7.753 | 5.37 |
| Supersmoother / `stats::supsmu` | 0.106 | 0.098 | 0.93 |
| Robinson + HC3 / R matrix implementation | 8.200 | 13.047 | 1.59 |

R's nonsmoothed `rq` took 0.387 ms and returned coefficient-sum checksum
1.848882736561; it is deliberately not the R denominator in the smoothed-loss
speed ratio. The matching smoothed fits gave 1.848624679258 (Java) and
1.848624649982 (R). Other checksums agreed to the printed precision:
OLS 2.034455174475; multinomial log likelihood -532.29220145647;
kernel 512.474170142173; supersmoother 512.985720853800;
Robinson 1.398956536347.

R's high-level model constructors may do more bookkeeping than the Java APIs;
`npreg` includes its normal result preparation. Robinson is compared to explicit
R matrix algebra, not to all features of `npplreg`. No end-to-end CSV/CLI startup,
bandwidth selection, or multivariate hypothesis testing is included. R's
supersmoother was about 7% faster on this fixture.

## Reproduce

From the repository root, after dependencies and benchmark classes are built:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/R/regression-accuracy.R
./gradlew.bat test --tests 'org.jlinalg.regression.*' --tests 'org.jlinalg.penalized.*'
./gradlew.bat '-Djlinalg.benchmark.regression.rows=512' `
  '-Djlinalg.benchmark.regression.warmups=100' `
  '-Djlinalg.benchmark.regression.measurements=101' benchmarkRegressionFamilies
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' `
  src/benchmark/resources/r-reference/benchmark-regression-families.R 512 7 100
```

R scripts prepend `build/r-library` and retain the user's existing libraries.
The Windows sandbox may require approval to read dependencies in the R user
library. The repaired benchmark consumes results; the pre-audit benchmark's
permanently zero checksum and false-convergence timings are superseded.
