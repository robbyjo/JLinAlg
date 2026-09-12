# Regression families

This vignette covers six regression families and their numerical contracts.
The examples use primitive arrays so the numerical contract is visible
without a dataframe dependency.

## Multivariate Gaussian regression

`MultivariateRegression` reuses one pivoted QR factorization across responses
and returns outcome-major coefficients, row-major fitted values and residuals,
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
class logits with column-scaled BFGS and a stable log-sum-exp likelihood:

```java
MultinomialRegressionResult fit = MultinomialRegression.fit(classes, design, 3);
double[] probabilities = fit.probabilities(); // row-major, class-major
```

The API is close to the common `nnet::multinom` contract, but it is not a
claim of optimizer-level equality. Inspect `logLikelihood()`, `iterations()`,
and `converged()` for every fit.

## Quantile regression

`QuantileRegression.fitExact` minimizes the actual nonsmoothed pinball loss
using a deterministic primal-dual linear program:

```java
QuantileRegressionResult fit = QuantileRegression.fitExact(y, design, 0.5);
double[] median = fit.coefficients();
QuantileLinearProgram.Result certified = QuantileLinearProgram.solve(y, design, 0.5);
double gap = certified.certificate().relativeGap();
```

The exact API checks primal/dual feasibility, complementarity, duality gap,
and original-unit reconstruction. `QuantileLinearProgram.Options` controls
iterations and tolerance. Full-column-rank designs with n >= p are required;
an intercept is an explicit design column. Tied/nonunique optima need not match
one arbitrary R coefficient vector. "Exact" identifies the nonsmoothed
objective, not exact arithmetic or quantile standard errors. Large offsets may
make the requested original-unit tolerance unrepresentable; inspect convergence.

The compatibility API `QuantileRegression.fit` still minimizes a smoothed
pinball loss, with `QuantileRegressionOptions` controlling smoothing and line
search. Its score-checked convergence and independent smoothed R comparison
remain separate from the exact LP. Neither API implements every nonparametric
or inferential feature of `quantreg`/`quantreg.nonpar`.
See [exact LP accuracy and paired R timing](../../src/benchmark/resources/exact-quantile-benchmark/exact-quantile-results.md).

### Covariance for exact quantile fits

```java
// f_i is the conditional response density at the fitted quantile, in 1/y units.
var inference = QuantileRegressionInference.fitExact(y, design, 0.5, densities);
double[] covariance = inference.coefficientCovariance(); // row-major
double[] standardErrors = inference.standardErrors();
// Alternative for iid errors only, with a caller-set response-unit bandwidth:
var iid = QuantileRegressionInference.fitExactIidKernel(y, design, 0.5, 0.4);
```

For A = X' diag(f_i) X the supplied-density covariance is
`tau*(1-tau) * A^-1 * X'X * A^-1`. Densities may vary across observations;
they must be positive, finite, and either known or consistently estimated.
Predictor and density scaling plus pivoted QR avoid an unscaled normal-equation
inverse; observation influence outer products form the covariance.
The iid alternative estimates the common residual density at zero with a pooled
Gaussian kernel. It is not a heteroskedastic conditional-density estimator.

Both paths require a converged exact fit, n > p, independent observations,
a correctly specified linear conditional quantile, and regular continuous
response densities. They do not provide valid mass-point or clustered-response
inference. These are asymptotic errors, not finite-sample exact inference or
Student-t errors. Bandwidth selection is intentionally separate: a positive
fixed bandwidth alone does not ensure valid inference; the usual shrinking
bandwidth and density regularity conditions are still needed. Singular designs,
failed fits, invalid densities, and numerically unresolved covariance are rejected.
See [R and analytic validation](../feasible-extensions-validation.md).

## Nonparametric regression and supersmoothing

`KernelRegression` provides one-dimensional Gaussian Nadaraya--Watson
regression. A bandwidth can be supplied explicitly or selected by its
Silverman-style default:

```java
KernelRegression.Result fit = KernelRegression.fit(x, y);
double[] prediction = fit.predict(new double[] {-1.0, 0.0, 1.0});
```

`SuperSmoother` implements Friedman's running-lines algorithm: smooth the
leave-one-out absolute errors for three spans, smooth the selected spans,
interpolate their fitted curves, and perform the final short-span smoothing:

```java
SuperSmoother.Result fit = SuperSmoother.fit(x, y);
double[] smooth = fit.fittedValues();
```

The algorithm is an attributed port of R's GPL-licensed `ppr.f` supersmoother.
Periodic fitting requires at least five observations; nonperiodic fitting requires four.
An overload exposes weights, fixed span (zero means CV), periodic fitting, and
bass control. Fitted values and selected spans retain input order; predictions
interpolate the fitted curve, with constant extrapolation at the boundaries.
Repeated queries are batch/order invariant. Missing and nonfinite input is
rejected rather than silently omitted; supplied weights must be positive.

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

When an intercept is present, the backfitting smooth is centered to mean zero.
For identifiable slopes with variation conditional on `z`, use Robinson
partialling-out with smoother-aware HC3-style slope covariance:

```java
PartiallyLinearInference.Result inference =
    PartiallyLinearRegression.fitWithInference(
        y, linearDesign, z, 0.4, BackendPolicy.CPU);
double[] slopeSe = inference.slopeStandardErrors();
```

The covariance concerns the reported `slopeIndices`, not the normalization
intercept. With smoother S, A=I-S, and residualized design Xt=AX, the coefficient
influence is L=(Xt'Xt)^-1 Xt'A and the residual map is B=A-Xt L. The sandwich
uses original-observation influence L and residual corrections e[i]/B[i,i];
it must not treat the smoothed residuals as independent OLS errors. This reduces
to ordinary HC3 when S=0, but is asymptotic and may be conservative for
high-leverage smoothing. It does not correct smoothing bias or bandwidth selection.

`homoskedasticSlopeCovariance()` alternatively returns (e'e / tr(BB')) LL'.
This is unbiased for homoskedastic noise when B annihilates the true mean;
otherwise smoothing bias remains. `residualNoiseDegreesOfFreedom()` returns
tr(BB'), not an exact Student-t denominator. A paired-observation noise test
checks the exact sampling variance and guards against the former factor-of-two
underestimate caused by ignoring the smoothing covariance.

## Accuracy and speed profile

`RegressionAccuracyTest` preserves the audit counterexamples.
`RegressionRParityTest` compares saved R fixtures for multivariate OLS,
`nnet::multinom`, smoothed quantile optimization and `quantreg::rq`, weighted,
fixed-span and periodic `stats::supsmu`, Gaussian kernel regression, Robinson
HC3, and `selectiveInference` truncation limits. Regenerate and test with:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/R/regression-accuracy.R
.\gradlew.bat test --tests 'org.jlinalg.regression.*'
.\gradlew.bat '-Djlinalg.benchmark.regression.rows=512' `
  '-Djlinalg.benchmark.regression.warmups=100' `
  '-Djlinalg.benchmark.regression.measurements=101' benchmarkRegressionFamilies
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' `
  src/benchmark/resources/r-reference/benchmark-regression-families.R 512 7 100
```

The benchmark consumes fitted results and rejects unconverged iterative fits.
R 4.6.1 is installed outside PATH on the validation host; reference scripts use
the project library at `build/r-library` plus installed user-library dependencies.
The two programs use identical deterministic inputs and emit comparable
checksums. Startup and reference-file I/O are excluded. The R quantile benchmark
separately reports the matching smoothed objective and the nonsmoothed `rq`
estimator: these must not be conflated. The semiparametric reference is explicit
R matrix algebra, not a claimed package-wide `npplreg` speed comparison.

## Boundaries

These APIs are deterministic and in-memory. They do not currently provide
formula parsing, missing-data policy, categorical encoding, bootstrap
intervals, or a general mixed-effect wrapper for the new six families. The
existing mixed, pedigree, mediation, beta, GLMM, and penalized APIs remain
the appropriate paths for those completed feature areas.
