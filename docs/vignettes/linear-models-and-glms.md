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
p-values treat the selected active set as fixed and ignore its selection.

### Selection-aware inference

`SelectionAwarePenalizedInference.fit` uses a deterministic held-out split:
only the first `selectionFraction` of rows select predictors; weighted OLS on
the remaining rows supplies inference. Full-data weights are split and
propagated. High-dimensional selection is allowed provided the selected
inference design has sufficient observations and rank. The rows in the two
parts must be independent; this is not a solution for clustered/time-dependent
rows split across the boundary. Do not select lambda using the inference rows.

For same-data Gaussian inference at an externally fixed penalty, use:

```java
PolyhedralSelectiveInference.Result inference =
    PolyhedralSelectiveInference.fit(y, predictors,
        0.12, 1.0,       // fixed lambda and alpha; alpha=1 is LASSO
        0.5, true,      // independently known noise SD; include intercept
        0.95, BackendPolicy.CPU);
```

The method derives the active-sign and inactive KKT inequalities, conditions
on that polyhedron, and inverts the truncated-normal distribution. It also
supports elastic net with `0 < alpha < 1`. Coefficient targets are the
selected-model least-squares projections, not shrunken penalized coefficients.
Predictors are not standardized; centering for the intercept is deterministic.
The selected design must have full column rank. Gaussian errors with known
independent noise SD and a response-independent lambda are essential: plugging
in a same-data residual SD or a CV-selected lambda does not preserve exactness.
Intervals are conditional on the selected signs as well as the active set.
Returned p-values are two-sided, whereas R's `fixedLassoInf` defaults to
sign-oriented one-sided tests. See Lee et al. (2016),
[Exact post-selection inference](https://arxiv.org/abs/1311.6238).

## Command-line regression

The same beta and penalized Gaussian fitters are available from the executable
JAR. Inputs are numeric CSV/TSV files with a header row. Both commands include
an intercept by default:

```console
java -jar jlinalg-<version>.jar beta-regression \
  --input proportions.tsv --response proportion --mean dose \
  --precision batch_score --out beta-coefficients.tsv

java -jar jlinalg-<version>.jar penalized-regression \
  --input continuous.tsv --response y --predictors x1,x2,x3 \
  --model elastic-net --alpha 0.5 --lambda-grid 1,0.3,0.1,0.03 \
  --cv-folds 5 --out penalized-coefficients.tsv
```

Use `--model ridge` or `--model lasso` for the endpoint penalties. Omit
`--precision` for constant beta-regression precision. The penalized command
reports the selected lambda, active-count, objective, and convergence flag;
`--no-intercept` and `--no-standardize` expose the corresponding API controls.
Cross-validation uses deterministic folds; pass `--seed` to change the fold
assignment.

## Formula equivalent

For an R-like user-facing layer, see the [formula vignette](formulas-and-backends.md).
Formula compilation happens once; numerical fitting still consumes contiguous
primitive arrays.
