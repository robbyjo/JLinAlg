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

Use the full Toeplitz Gaussian covariance for stationary series, missing
observations, or independent panels sharing parameters:

```java
ExactArmaResult exact = ExactArma.fit(
    seriesWithNaN, ArimaOrder.arma(1, 1), true,
    BackendPolicy.PREFERRED);

ExactArmaResult panel = ExactArma.fitPanel(
    List.of(siteOne, siteTwo, siteThree),
    ArimaOrder.ar(1), true, BackendPolicy.PREFERRED);
```

Missing values are exactly marginalized for stationary ARMA. Integrated models
are not accepted by `ExactArma`; JLinAlg does not currently claim a diffuse
Kalman likelihood for them. Complete stationary series use a Durbin-Levinson
likelihood; missing series retain the dense marginal-covariance path.

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
