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
minimum-norm coefficient as unique. Non-estimable coordinates have `NaN`
standard errors/tests/intervals, and `testContrast()` rejects contrasts outside
the design row space. The minimum-norm covariance remains available for
estimable combinations. Full-rank QR equilibrates extreme or disparate column
units, including mixed `[1e-6, 1e6]` units. Deficient fits identify rank in scaled
coordinates but retain the original-coordinate minimum norm. Their original-unit
SVD must reproduce the scaled reference's fitted projection and leverage;
inaccurate or nonconvergent factorizations raise an explicit reparameterization
error instead of silently dropping an estimable intercept. Very ill-conditioned
and some wide deficient designs can therefore be rejected rather than fitted.

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
quasi-Poisson families. Fixed-dispersion GLM coefficient tests use Wald z tests;
estimated dispersion uses residual Student t inference and F contrasts.
`associationStatistics()` follows the same distribution as `pValues()`.
Gaussian identity fitting uses one weighted least-squares solve, including
weights and offsets, rather than redundant IRLS factorizations.

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

Pearson dispersion scales the coefficient covariance, but it is not generally
the dispersion used for the likelihood. Built-in Gaussian and inverse-Gaussian
likelihoods use ML dispersion; Gamma numerically profiles its actual density.
The latter is not the approximate Gamma AIC calculation in R's `stats::glm`.
Custom families retain their supplied-dispersion likelihood contract. AIC
counts numerical rank, plus dispersion only when the family estimates it.
An iteration limit or a small accepted line-search step without stationary
final scoring equations is not reported as convergence. Binomial means below
`1e-12` and small positive Gamma/inverse-Gaussian responses are no longer
flattened by arbitrary response-scale variance floors. Floating-point overflow,
unrepresentable probabilities, separation, and near-singular models still
require care; a converged flag is not an identification certificate.
Binomial score, Fisher information, deviance and density use the predictor's
small tail directly, even when the reported mean rounds to one. Compatible
`GlmFamily` predictor-aware defaults let custom families preserve their existing
contracts. Information is not floored to fabricate finite standard errors:
unrepresentable working precision is an error. Exactly zero-residual Gaussian
ML has likelihood supremum `+Infinity` and AIC `-Infinity`.

## Clustered marginal models (GEE)

`Gee.fit()` supports weights, link-scale offsets, cluster IDs and repeated-wave
indices. Select the working correlation with `GeeOptions`; select model-based,
cluster-sandwich or a documented finite-cluster correction separately from the
working correlation. The ordinary sandwich does not imply small-sample exactness.
Use `BackendPolicy.CPU` for reproducible portable CPU timings; strict native/GPU
policies retain their existing availability checks. No native CHOLMOD runtime
certification follows from these portable CPU checks.

Changing the measurement units of a Gamma response must shift the log-link
intercept, not change slopes, dispersion or sandwich covariance. Regression
tests now cover this for independence and exchangeable working correlations.
Gaussian dispersion and naive covariance now scale with the response units,
including responses around `1e-12`; zero-residual fits report zero dispersion
and covariance. The score convergence norm is divided by the square root of
each sensitivity diagonal so tiny physical units do not prevent convergence.
Serial score accumulation reduces cluster contributions immediately instead of
retaining every cluster's bread/meat arrays. Parallel and serial results are
checked independently. The R audit uses an independently constructed exact GLS
solution and cluster sandwich for fixed correlation; it is not a `geepack`
benchmark or validation of every estimated-association/bias-adjustment variant.

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

Observation weights are normalized without summing huge unscaled values.
Standardization works for tiny nonconstant predictor units; it does not drop a
column merely because its SD is below `1e-14`. Near-exact fits reconstruct RSS
from residuals when the covariance-form quadratic loses precision.

CV now requires convergence of every full-data and training-fold fit. Its mean
risk is the held-out weighted MSE across observations, not an unweighted average
of unequal folds. The grouped SE is the square root of the fold-weighted mean
squared deviation divided by `K-1`; folds remain deterministic for a given seed.
Automatic CV reuses its already fitted full-data path. Tune the iteration budget
and tolerance if a path is rejected; nonconverged fits are not eligible winners.
CV preserves raw weights for each training subset. Evaluation weights smaller
than representable relative mass contribute zero, without triggering a false
positivity error or contaminating risk with `0/0` from a zero-mass fold.
Penalty factors in this API multiply the stated L1/L2 terms directly and are
not silently rescaled to sum to the predictor count as in `glmnet`.

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

## Reproducing the v0.3.0 fitting audit

The independent generator is `src/benchmark/r/fitting_audit_v030.R`; frozen
fixtures and seed-17 fold assignments are under
`src/test/resources/fitting-audit-v030`. It uses R 4.6.1 and `glmnet` 5.0, base R
weighted/offset GLMs, independently profiled Gamma/inverse-Gaussian densities,
and exact active-sign enumeration for small ridge/elastic-net fits with an
unpenalized predictor. It installs nothing. On Windows, reading installed R
user-library dependencies may require an elevated reference process.

`org.jlinalg.benchmark.FittingAuditBenchmark` and that R script construct the same
3,000-row data. Each timing is a median of seven batches of 20 warmed fits, with
convergence checks and consumed coefficient checksums. Measured portable-Java
CPU times were 0.523 ms OLS, 2.103 ms Poisson, 4.577 ms fixed-correlation GEE,
and 0.204 ms for a three-lambda weighted LASSO path. Corresponding R medians
were 0.5, 6, 22.5, and 1 ms; R timings have 0.5 ms batch resolution on this host.
OLS/GLM Java returns full inference while the timed R fits do not materialize
every summary; Java GEE also computes auxiliary diagnostics. These are narrow
workload timings, not interchangeable feature-cost or universal speed claims.
Checksums agree to floating-point rounding. Full evidence, allocation counts
and limitations are in `src/benchmark/resources/fitting-audit-v030/evidence.md`.
