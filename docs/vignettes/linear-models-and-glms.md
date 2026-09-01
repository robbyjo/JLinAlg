# Linear models, GLMs, and penalized regression

## Ordinary least squares

Include an intercept explicitly as a column of ones:

```java
double[] y = {1, 2, 5, 7, 9};
double[][] x = {
    {1, 0},
    {1, 1},
    {1, 2},
    {1, 3},
    {1, 4}
};

OlsResult fit = Ols.fit(y, x);
System.out.println(Arrays.toString(fit.coefficients()));
System.out.println(Arrays.toString(fit.standardErrors()));
System.out.println(Arrays.toString(fit.tStatistics()));
System.out.println(Arrays.toString(fit.pValues()));
```

`fit.coefficients()[1]` is the change in the response per unit of the second
design column. OLS uses `N-rank(X)` residual degrees of freedom. Inspect
`fit.rank()` and the backend provenance when results will be persisted.

For a rank-deficient model, the default is a clear error. A minimum-norm SVD
solution is an explicit opt-in:

```java
OlsOptions options = new OlsOptions(
    RankDeficiencyStrategy.MINIMUM_NORM, 0.95,
    MissingDataPolicy.ERROR);
OlsResult fit = Ols.fit(y, x, options, BackendPolicy.PREFERRED);
```

The resulting individual coefficients depend on the identifying convention;
test scientifically meaningful estimable contrasts rather than treating every
minimum-norm coefficient as unique.

## Weights, offsets, and missing rows

Prior weights are positive and offsets are additive on the response scale for
OLS. `MissingDataPolicy.OMIT` compacts complete rows once and preserves their
original indices in the result:

```java
double[] weights = {1, 1, 2, 2, 1};
double[] offset = {0, 0, 0.2, 0.2, 0.2};
OlsOptions options = new OlsOptions(
    RankDeficiencyStrategy.ERROR, 0.95,
    MissingDataPolicy.OMIT);

OlsResult weighted = Ols.fit(
    y, x, weights, offset, options, BackendPolicy.PREFERRED);
int[] retainedRows = weighted.retainedRows();
```

Use `testContrast()` for a joint hypothesis. The rows of the contrast matrix
select linear combinations of coefficients.

## Generalized linear models

The family defines the variance and link; the design convention is unchanged:

```java
double[] disease = {0, 0, 0, 1, 1, 1};
double[][] design = {
    {1, -2}, {1, -1}, {1, 0},
    {1, 0.5}, {1, 1}, {1, 2}
};

GlmResult logistic = Glm.fit(
    disease, design, GlmFamilies.binomial());
System.out.println(logistic.coefficients()[1]);
System.out.println(logistic.pValues()[1]);
```

Available factories include Gaussian, binomial, Poisson, Gamma,
inverse-Gaussian, fixed-size negative-binomial, quasi-binomial, and
quasi-Poisson families. GLM coefficient tests are asymptotic Wald z tests.

For grouped binomial data, supply response proportions and trial counts as
prior weights. For count rates, use log exposure as an offset:

```java
double[] counts = {3, 7, 4, 12};
double[][] countDesign = {{1, 0}, {1, 1}, {1, 0}, {1, 1}};
double[] exposure = {10, 20, 10, 30};
double[] logExposure = Arrays.stream(exposure).map(Math::log).toArray();

GlmResult poisson = Glm.fit(counts, countDesign,
    GlmFamilies.poisson(), null, logExposure,
    GlmOptions.defaults(), BackendPolicy.PREFERRED);
```

Check `converged()`, deviance, Pearson residuals, dispersion, and AIC. AIC is
not defined for quasi-likelihood families.

## Ridge, LASSO, and elastic net

Penalized predictors should normally exclude the intercept; the fitter handles
its unpenalized intercept separately:

```java
double[][] predictors = {
    {0, 1}, {1, 0}, {2, 1}, {3, 0}, {4, 1}
};

PenalizedRegressionResult ridge =
    PenalizedRegression.ridge(y, predictors, 0.1);
PenalizedRegressionResult lasso =
    PenalizedRegression.lasso(y, predictors, 0.1);

ElasticNetOptions enet = ElasticNetOptions.builder()
    .alpha(0.5)                 // 0=ridge, 1=LASSO
    .standardize(true)
    .build();
PenalizedRegressionResult elastic =
    PenalizedRegression.fit(y, predictors, 0.1, enet);
```

Choose lambda with reproducible K-fold cross-validation rather than selecting
it on the final analysis response:

```java
PenalizedCrossValidationResult cv =
    PenalizedRegressionCrossValidation.automatic(
        y, predictors, 50, 1e-3, 5, 12345L, enet);
double lambdaMin = cv.lambdaMinimum();
double lambdaOneSe = cv.lambdaOneStandardError();
```

Use `PenalizedRegressionInference.ridge` for model-based ridge inference.
For LASSO/elastic net, `refitActiveSet` performs an optional OLS refit, but its
p-values are conditional on the selected active set and are not automatically
selection-valid.

## Formula equivalent

For an R-like user-facing layer, see the [formula vignette](formulas-and-backends.md).
Formula compilation happens once; numerical fitting still consumes contiguous
primitive arrays.
