# Time-series models and ARIMA-error mixed models

## AR, MA, and ARMA

```java
ArimaResult ar1 = Arima.fit(series, ArimaOrder.ar(1));
ArimaResult ar3 = Arima.fit(series, ArimaOrder.ar(3));
ArimaResult ma1 = Arima.fit(series, ArimaOrder.ma(1));
ArimaResult arma = Arima.fit(series, ArimaOrder.arma(2, 1));

if (!arma.converged()) {
    throw new IllegalStateException(arma.convergenceMessage());
}
System.out.println(Arrays.toString(arma.autoregressive()));
System.out.println(Arrays.toString(arma.movingAverage()));
System.out.println(arma.aicc());
```

MA signs match R: a positive MA coefficient enters the observation equation
with a positive sign. The conditional fitter transforms parameters to preserve
stationarity and invertibility. Pure nonseasonal AR models use their
closed-form conditional least-squares solution. Other models use one
deterministic optimizer start by default; request additional starts for a
particularly irregular likelihood:

```java
ArimaOptions robust = ArimaOptions.builder()
    .optimizationStarts(5)
    .build();
```

## Integrated and seasonal models

```java
ArimaOptions seasonal = ArimaOptions.builder()
    .seasonalOrder(SeasonalArimaOrder.of(1, 1, 1, 12))
    .build();

ArimaResult sarima = Arima.fit(
    monthlySeries, ArimaOrder.arima(1, 1, 1), seasonal);
ArimaForecast nextYear = sarima.forecast(12, 0.95);
```

For exactly one ordinary or seasonal difference, represent deterministic drift
explicitly:

```java
ArimaOptions drift = ArimaOptions.builder()
    .includeDrift(true)
    .build();
ArimaResult randomWalk = Arima.fit(
    series, ArimaOrder.arima(0, 1, 0), drift);
```

`includeMean` applies only to undifferenced models. Forecast means are returned
on the original scale with innovation-based normal intervals.

## Diagnostics

```java
double[] acf = arma.residualAutocorrelation(20);
LjungBoxResult lb = arma.ljungBox(20);
double[] pacf = TimeSeriesDiagnostics.partialAutocorrelation(
    arma.innovations(), 20);

System.out.printf("Ljung-Box=%g df=%d p=%g%n",
    lb.statistic(), lb.degreesOfFreedom(), lb.pValue());
```

Choose diagnostic lags before inspecting the result and account for the fitted
ARMA parameter count. Residual plots and domain-specific intervention checks
remain the caller's responsibility.

## Automatic small-order search

```java
ArimaSelectionResult selected = AutomaticArima.select(
    series, 3, 1, 3);
ArimaResult best = selected.bestModel();
selected.candidates().forEach(candidate ->
    System.out.println(candidate.order() + " " + candidate.aicc()));
```

This is exhaustive AICc selection over the requested nonseasonal bounds. It is
not a stepwise forecasting oracle; review convergence, residual diagnostics,
and scientifically plausible orders.

## Exact stationary ARMA

Use the full stationary Gaussian likelihood for complete series, missing
observations, or independent panels sharing parameters:

```java
ExactArmaResult exact = ExactArma.fit(
    seriesWithNaN, ArimaOrder.arma(1, 1), true,
    BackendPolicy.PREFERRED);

ExactArmaResult panel = ExactArma.fitPanel(
    List.of(siteOne, siteTwo, siteThree),
    ArimaOrder.ar(1), true, BackendPolicy.PREFERRED);
```

Missing values are exactly marginalized for stationary ARMA. Both complete and
missing series use state-space innovations with an unconditional stationary
initial covariance. This has the same Gaussian likelihood, innovation-variance
scale, mean interpretation, AIC/BIC conventions and coefficient covariance API
as the previous stationary implementation. Independent panels reset the state
for each series and share the estimated parameters and innovation variance.
`ExactArma` continues to reject integrated orders.

Only NaN denotes missingness; infinite observations are rejected. Positive
innovation variances are accepted without an absolute scale cutoff. Degenerate
zero-variance likelihoods raise an exception rather than returning an optimizer
penalty as a fitted likelihood. Conditional AR estimates are optional starting
values, so a short series valid for exact fitting is not rejected merely because
conditional fitting requires more observations.

Check `exact.coefficientInferenceAvailable()` before using coefficient SEs.
Observed information must be positive definite beyond finite-difference error
at two step sizes. Singular or unresolved information produces NaN covariance
and SEs; no ridge is added to make an unidentified AR/MA decomposition appear
estimable. `converged()` describes optimization and is separate from this
inference diagnostic. Nonconverged estimates also have unavailable covariance.

`SparseMissingSeries.fit` uses this same state-space implementation. Its path
metadata describes state-space innovations; it no longer claims an observed-
pattern sparse Cholesky factorization. No original-series or observed-series
quadratic covariance array is allocated. The filter uses CPU arithmetic;
coefficient-Hessian linear algebra follows the requested backend policy.

## Exact diffuse integrated likelihood

```java
DiffuseArima.Result exactIntegrated = DiffuseArima.fit(
    seriesWithNaN, ArimaOrder.arima(0, 1, 1),
    ArimaOptions.builder()
        .seasonalOrder(SeasonalArimaOrder.of(0, 1, 1, 12))
        .optimizationStarts(3)
        .build());
ArimaResult fitted = exactIntegrated.fit();
if (!fitted.converged()) {
    throw new IllegalStateException(fitted.convergenceMessage());
}
ArimaForecast forecast = fitted.forecast(12);
```

`DiffuseArima` now fits the **original observation sequence** by exact diffuse
Kalman recursions. It supports ordinary integration, seasonal-only integration,
their combination, multiplicative seasonal AR/MA, and arbitrary NaN observation
patterns. Infinity is rejected. `Arima.fit` remains the separate conditional
fitter shown earlier; selecting `DiffuseArima` does not call that fitter.

Write the differencing polynomial as
`delta(B) = (1-B)^d (1-B^s)^D`, with degree `b = d+sD`. The state consists of a
stationary innovations-form ARMA block of size `r = max(p_effective,q_effective+1)`
and `b` integrated lags. The observation vector is
`Z = [1, 0, ..., 0, -delta_1, ..., -delta_b]`. The stationary transition has AR
coefficients in its first column and ones on its superdiagonal. The first
integrated lag receives `Z`, and remaining integrated lags shift. The process
noise loading is `[1, theta_1, ..., theta_q, 0, ...]` in the stationary block.

The initial covariance is represented symbolically as `P_* + kappa P_inf`,
where `P_inf` is identity on the integrated lags and zero elsewhere. The
stationary finite block solves `P = T P T' + R R'` by matrix doubling; it does
not truncate a fixed-length impulse response. The filter starts at the first
observed position. An entirely missing prefix
is integrated out analytically by shifting the diffuse time origin; its length
cannot accumulate cancelling covariance terms or consume diffuse rank. The
original time grid is retained for residuals, trends and forecasts. At each
observed update, let `u_* = P_* Z'`, `u_inf = P_inf Z'`, `F_* = Z u_*`,
and `F_inf = Z u_inf`.
When `F_inf > 0`, the exact measurement updates are:

```text
a      <- a + u_inf * v / F_inf
P_inf  <- P_inf - u_inf u_inf' / F_inf
P_*    <- P_* - (u_inf u_*' + u_* u_inf') / F_inf
                 + F_* u_inf u_inf' / F_inf^2
```

When `F_inf = 0`, the ordinary finite Kalman update applies. Missing positions
receive prediction only; they do not decrement the diffuse rank. Fitting
rejects a pattern that leaves diffuse directions unidentified (for example,
never observing one season of a seasonal random walk). `diffuseStateCount()`
is `d+sD`, not the number of differencing operators. `diffuseLikelihood()`
identifies this algorithm and is separate from optimizer convergence.

Likelihood normalization follows [R's `stats::arima(method="ML")`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/arima.html):
diffuse observations contribute neither squared residuals nor determinants.
With `m` finite innovations, `Q = sum(v^2/F_*)` and `L = sum(log(F_*))`, the
reported innovation variance is `Q/m` and
`logLik = -0.5 * [m * (log(2*pi) + 1 + log(Q/m)) + L]`.
The mathematical diffuse-density normalization can additionally include
`-0.5 sum(log(F_inf))`; that term is independent of the ARMA parameters for this
initialization and observation pattern and is omitted here to retain R's
convention. Likelihoods from different missing patterns or differencing orders
should not be compared without considering this normalization.

R uses finite `kappa=1e6` initialization and a diffuse-gain cutoff, whereas Java
tracks the infinite and finite coefficients separately. Consequently agreement
with R is approximate for integrated models. The frozen reference files also
contain likelihoods with several values of `kappa`. At `kappa=1e7`, all 11
integrated fixtures agree with Java within `8.7e-6` at the same coefficients;
the largest discrepancy at R's default `kappa=1e6` is `8.7e-5`. Small `kappa`
can change which observations R regards as diffuse and is not a reliable
reference for the limiting likelihood.

Forecasts propagate the final filtered state and covariance, including
uncertainty from missing terminal observations and stationary initial states.
They do not substitute missing values or feed conditional residuals into a
forecast recursion. Intervals are conditional on fitted coefficients and
innovation variance; parameter-estimation uncertainty is not included.
For a single differencing operator, `includeDrift(true)` estimates the mean of
the differenced process. With a seasonal difference of period `s`, this is `s`
times the level-scale linear trend slope. `includeMean` is ignored for
integrated models. Drift with multiple differencing operators is rejected.

Drift is profiled by unrestricted GLS through the Kalman innovations at each
dynamic-parameter candidate. The observation and time-trend series share the
same gains; recursive least squares accumulates the profiled residual sum of
squares without subtracting two nearly equal quadratic forms. No bounds based
on the small subset of observable explicit differences can exclude the drift
MLE. Profiling removes drift from the nonlinear search but retains its parameter
count in AIC, AICc and BIC.

One-dimensional optimization samples 64 bounded intervals plus the supplied
initial point, then refines each sampled local minimum. It retains the best
seed/boundary value and checks the final projected score; it does not assume
that the entire ARMA likelihood is unimodal. Evaluation limits are respected.
This is a multi-basin numerical search, not a proof of global optimality for
arbitrarily narrow unsampled minima. Higher-dimensional fits can still benefit
from multiple starts.

Exact innovations are `v/sqrt(F_*)` on the original time grid, with NaN for
missing and diffuse observations. ACF and Ljung-Box methods reject internal
missing/diffuse gaps rather than silently collapse time. `differencedSeries()`
is available for inspection but is not used to fit or forecast the exact model.

State storage is `O((r+b)^2)`, plus `O(n)` input/result arrays. Sparse transition
rows and symmetric covariance updates give `O(n (r+b)^2)` filtering work.
Stationary initialization uses state-sized matrix doubling. High orders and
long seasonal periods can still be expensive or ill-conditioned; this is a
covariance filter, not a square-root filter. The diffuse pivot tolerance is
`1e-9` in unit diffuse-state coordinates. Unrestricted regression terms,
coefficient covariance/standard errors for `DiffuseArima`, and smoothing of
missing historical values are not exposed by this API. The stationary
`ExactArmaResult` covariance/SE API remains available.

## Accuracy and matching-method timings

The frozen R fixtures and generator are in
`src/test/resources/timeseries/`. They cover complete/missing stationary ARMA,
ARIMA(1,1,0), ARIMA(0,1,1), missing ARIMA(1,1,1), seasonal-only integration,
ordinary plus seasonal differencing, twice-integrated MA, ordinary and seasonal
drift, multiplicative ordinary/seasonal ARMA, and monthly missing observations.
Each contains coefficients, innovation variance, ML likelihood, effective
sample size, and eight-step forecast means and standard errors. Generation
uses R 4.6.1, base `stats`, seeded simulations, `method="ML"`,
`SSinit="Rossignol2011"`, and `reltol=1e-12`; no optional package is needed.

`ExactDiffuseArimaTest` checks both fixed-parameter filtering and optimized
fits. Independent tests cover analytical random-walk gaps, complete-series
equivalence to exact likelihood after differencing, a dense stationary
marginalization oracle, a stationary AR root of `0.999999`, unidentified
seasonal directions, budget-exhausted optimizers, repeatable forecasts, and a
50,000-position missing series. The independent-review regression suite also
covers 100,000 leading NaNs, ordinary and seasonal missing-data drift MLEs,
competing MA(1) likelihood basins, scales from 1e-50 to 1e50, a short exact AR(2),
unidentified AR/MA cancellation, and both signs of infinity. Its additional R
references are generated by `generate-review-arima.R` in the same fixture
directory. All 35 time-series package tests passed in an
isolated JVM with `-Xmx128m` on 2026-09-08.

Two paired benchmark sources distinguish the work being compared:

- `ExactArimaLikelihoodBenchmark` and `exact_arima_likelihood_benchmark.R`
  time fixed-parameter likelihood evaluation **including stationary state
  initialization**, excluding optimization, coefficient Hessians, forecasts,
  and input loading. Both consume every likelihood. The R benchmark calls the
  installed `stats` state constructor and Kalman likelihood at the same frozen
  coefficients. It uses R's finite-diffuse approximation.
- `ExactTimeSeriesBenchmark` and `exact_time_series_benchmark.R` time converged
  full ML fits and eight-step forecasts for integrated cases, with accuracy
  gates and consumed checksums. Their optimizer work differs: R computes a
  coefficient covariance, while the diffuse Java API currently does not.
  Stationary Java fits use five starts; integrated Java fits use one. These are
  end-to-end timings, not an optimizer or inference-work equivalence claim.

Measured results and exact reproduction commands are recorded in
`src/benchmark/resources/timeseries/exact-arima-measurements.md`. Measurements
are local wall-clock medians and are not portable speed guarantees.

## LMM with ARIMA errors

```java
ArimaErrorLmmResult correlated = ArimaErrorLinearMixedModel.fit(
    response, fixedDesign, randomTerms,
    ArimaOrder.arma(1, 1),
    ArimaErrorLmmOptions.builder()
        .remlOptions(RemlOptions.defaults())
        .build(),
    BackendPolicy.PREFERRED);
```

The fitter profiles mixed-model variances by REML at each stationary error-
correlation candidate. When differencing is requested, it differences the
response and every fixed/random design column together. A level intercept then
becomes zero and is rejected; use a level-scale time trend to represent drift.

The current API treats rows as one ordered series. Use `ExactArma.fitPanel` for
block-independent shared ARMA parameters outside an LMM.
