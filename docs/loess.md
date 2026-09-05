# LOESS local polynomial regression

`Loess` implements one-dimensional direct-surface local polynomial regression.
It is deliberately separate from the spline/penalty GAM APIs: LOESS selects a
nearest-neighbor window for every evaluation point and fits a tricube-weighted
polynomial centered at that point.

```java
LoessOptions options = new LoessOptions(
    0.6,                  // span: fraction of observations in each neighborhood
    2,                    // local polynomial degree: 0, 1, or 2
    LoessFamily.GAUSSIAN,
    4);                   // R-compatible robust-pass convention

LoessResult fit = Loess.fit(x, y, priorWeights, options);
double[] fitted = fit.fittedValues();
double[] residuals = fit.residuals();
double trace = fit.effectiveDegreesOfFreedom();
double[] predictions = fit.predict(new double[] {-1.0, 0.0, 1.0});
```

Inputs must be finite, `0 < span <= 1`, and the predictor must vary. Prior
weights may be zero but not negative. The direct predictor extrapolates by a
local polynomial outside the observed range; such extrapolation can be
unstable and should be interpreted cautiously.

## Robust symmetric fitting

`LoessFamily.SYMMETRIC` applies Tukey-bisquare residual weights. The default
four-pass convention matches `stats::loess(..., family="symmetric",
control=loess.control(iterations=4))`. The final weights are available through
`robustnessWeights()` so complete rejection of influential points is visible.

For repeated response fits at one predictor geometry, sort once:

```java
Loess.Prepared prepared = Loess.prepare(x, options);
LoessResult first = prepared.fit(y1, weights);
LoessResult second = prepared.fit(y2, weights);
```

## R accuracy contract

The regression fixture uses R 4.6.1 `stats::loess` with
`surface="direct"`, `statistics="exact"`, and `trace.hat="exact"`.
Checked values cover weighted degree-two fitted values, arbitrary predictions,
trace of the smoother matrix, robust fitted values, and robustness weights.
The committed tolerances are `2e-12` for Gaussian output and `2e-10` for the
iteratively robust output.
Regenerate the source values with
`src/test/resources/r-reference/generate-loess-reference.R`.

This is not R's interpolated kd-tree surface. Current scope is one numeric
predictor, degree zero through two, finite complete data, direct prediction,
and fitted-point diagonal leverage. Multivariate neighborhoods, standard-error
prediction, NA exclusion/reinsertion, parametric columns, and R's interpolated
surface are not silently approximated.

## Measured performance

The deterministic 2026-09-05 benchmark used 5,000 evenly spaced rows, span
0.2, degree two, Gaussian fitting, two warm-ups, and five measurements on the
Windows 11 Intel Core i9-9900K development host:

| Runtime | Median seconds | Rows/second | Speedup vs R |
| --- | ---: | ---: | ---: |
| R `stats::loess`, direct | 0.090000 | 55,556 | 1.00x |
| JLinAlg end-to-end | 0.023688 | 211,081 | 3.80x |
| JLinAlg prepared geometry | 0.022872 | 218,606 | 3.93x |

R was configured with `statistics="none"` for the timed fit. JLinAlg still
computed fitted values, residuals, diagonal leverage, and exact trace, making
the timing conservative for JLinAlg. Both checksums were
`0.001584182714`.

Reproduce with:

```powershell
.\gradlew.bat benchmarkLoess
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\loess_benchmark.R 5000 0.2 5 2
```
