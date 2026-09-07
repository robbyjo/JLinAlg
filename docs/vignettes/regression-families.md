# Regression families

This vignette covers the six analysis families completed in the current TODO
review. The examples use primitive arrays so the numerical contract is visible
without a dataframe dependency.

## Multivariate Gaussian regression

`MultivariateRegression` fits one OLS system per response using a shared
predictor matrix and returns row-major coefficients, fitted values, residuals,
and the residual covariance matrix:

```java
MultivariateRegressionResult fit = MultivariateRegression.fit(response, design);
double[] beta = fit.coefficients(); // outcome-major, then predictor
double[] sigma = fit.residualCovariance();
```

This is the direct Gaussian reference case for `lm`/`rrr`-style shared-design
workflows. It requires a full-column-rank design and uses the configured JLinAlg
compute backend.

## Multinomial logistic regression

`MultinomialRegression` uses class `0` as the baseline and estimates the other
class logits with a deterministic, monotone line-searched gradient ascent:

```java
MultinomialRegressionResult fit = MultinomialRegression.fit(classes, design, 3);
double[] probabilities = fit.probabilities(); // row-major, class-major
```

The API is close to the common `nnet::multinom` contract, but it is not a
claim of optimizer-level equality. Inspect `logLikelihood()`, `iterations()`,
and `converged()` for every fit.

## Quantile regression

`QuantileRegression` minimizes a smoothed pinball loss, which keeps the fit
deterministic and avoids bootstrap sampling:

```java
QuantileRegressionResult fit = QuantileRegression.fit(y, design, 0.5);
double[] median = fit.coefficients();
```

The smoothing parameter and line-search controls are exposed through
`QuantileRegressionOptions`. This is a stable first implementation; it is not
an exact reproduction of every simplex, Frisch--Newton, or nonparametric path
in `quantreg` and `quantreg.nonpar`.

## Nonparametric regression and supersmoothing

`KernelRegression` provides one-dimensional Gaussian Nadaraya--Watson
regression. A bandwidth can be supplied explicitly or selected by its
Silverman-style default:

```java
KernelRegression.Result fit = KernelRegression.fit(x, y);
double[] prediction = fit.predict(new double[] {-1.0, 0.0, 1.0});
```

`SuperSmoother` selects among candidate spans independently at each training
point with leave-one-out local-linear tricube fits:

```java
SuperSmoother.Result fit = SuperSmoother.fit(x, y);
double[] smooth = fit.fittedValues();
```

The span selector is inspired by `stats::supsmu`; it deliberately exposes the
selected spans because edge behavior and candidate-span choices are part of
the result. It is not presented as a hidden-compatible clone of R's Fortran
implementation.

## Semiparametric partially linear regression

`PartiallyLinearRegression` estimates

\[
y_i = x_i^T\beta + g(z_i) + \epsilon_i
\]

by alternating least-squares updates for `beta` and Gaussian-kernel updates
for `g`:

```java
PartiallyLinearRegression.Result fit =
    PartiallyLinearRegression.fit(y, linearDesign, nonlinearPredictor, 0.5);
double[] beta = fit.coefficients();
double[] g = fit.smoothEffect();
```

This is a useful bounded semiparametric contract, not a package-wide claim
about all spline, single-index, varying-coefficient, or generalized additive
semiparametric models. The returned convergence flag and residuals should be
recorded with the analysis.

## Accuracy and speed profile

`RegressionFamiliesTest` checks coefficient recovery for multivariate OLS,
probability normalization for multinomial regression, and finite output for
the remaining four families. Run:

```powershell
.\gradlew.bat test --tests org.jlinalg.regression.RegressionFamiliesTest
.\gradlew.bat benchmarkRegressionFamilies `
  -Pjlinalg.benchmark.regression.rows=512 `
  -Pjlinalg.benchmark.regression.measurements=5
```

The benchmark prints median seconds and rows/second for each Java estimator.
Direct R timing was not available on the development host because `Rscript`
was not installed, so this release does not claim an R speedup for these six
new families. When R is available, compare the Gaussian multivariate case to
`lm`/`rrr`, multinomial fits to `nnet::multinom`, supersmoothing to
`stats::supsmu`, quantiles to `quantreg`, and kernel fits to the selected
`np`/`npreg` contract using identical data, tolerances, and warm-up policy.

## Boundaries

These APIs are deterministic and in-memory. They do not currently provide
formula parsing, missing-data policy, categorical encoding, bootstrap
intervals, or a general mixed-effect wrapper for the new six families. The
existing mixed, pedigree, mediation, beta, GLMM, and penalized APIs remain
the appropriate paths for those completed feature areas.
