# Individual-level instrumental-variable regression

Use two-stage least squares (2SLS) when a linear outcome model contains a
predictor that is correlated with the structural error, and you have excluded
instruments that shift that predictor without directly shifting the outcome.
This is an individual-level estimator; summary-association Mendelian
randomization uses the separate `org.jlinalg.mr` APIs.

An IV estimate is causal only under substantive assumptions that software
cannot verify:

1. the excluded instruments predict the endogenous regressors (relevance);
2. the instruments are independent of unmeasured outcome causes, conditional
   on the included exogenous regressors (exchangeability);
3. the instruments affect the outcome only through the modeled endogenous
   regressors (exclusion); and
4. the linear structural model and observation/cluster independence assumptions
   used for inference are appropriate.

## Fit a model

The example below estimates the effect of one endogenous exposure. The
exogenous matrix contains the intercept and an observed covariate. The excluded
instrument matrix contains two instruments. Include the intercept explicitly
when the model needs one.

```java
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.regression.InstrumentalVariableCovariance;
import org.jlinalg.regression.InstrumentalVariableOptions;
import org.jlinalg.regression.InstrumentalVariableRegression;
import org.jlinalg.regression.InstrumentalVariableResult;

double[] outcome = {2.2, 3.1, 2.7, 5.8, 4.9, 6.4, 5.7, 7.6};
double[][] exogenous = {
    {1, 0}, {1, 1}, {1, 0}, {1, 1},
    {1, 0}, {1, 1}, {1, 0}, {1, 1}
};
double[][] endogenous = {
    {0.7}, {1.2}, {1.0}, {2.1}, {1.8}, {2.6}, {2.4}, {3.2}
};
double[][] excludedInstruments = {
    {-1.2,  0.3}, {-0.7, -0.8}, {-0.3,  1.0}, {0.2, -0.4},
    { 0.5,  0.7}, { 0.9, -1.1}, { 1.3,  0.4}, {1.8, -0.2}
};

InstrumentalVariableOptions options = new InstrumentalVariableOptions(
    InstrumentalVariableCovariance.HC1, 0.95);
InstrumentalVariableResult fit = InstrumentalVariableRegression.fit(
    outcome, exogenous, endogenous, excludedInstruments,
    options, BackendPolicy.CPU);

// Coefficient order is [all exogenous columns, all endogenous columns].
double exposureEffect = fit.coefficients()[2];
double exposureSe = fit.standardErrors()[2];
double exposureP = fit.associationStatistics().pValues()[2];
double lower = fit.confidenceLower()[2];
double upper = fit.confidenceUpper()[2];
```

`HC1` is heteroskedasticity-robust and uses asymptotic normal inference. Choose
`HOMOSKEDASTIC` only when a constant structural-error variance is defensible;
it uses residual Student-t inference. `HC0` returns the unadjusted observation
sandwich.

The coefficient fit is

\[
\hat\beta=(X'P_ZX)^{-1}X'P_Zy,
\]

where `X` contains included exogenous and endogenous regressors, while `Z`
contains included exogenous regressors and excluded instruments. JLinAlg
computes the QR projections without forming the dense `P_Z` matrix. Fitted
outcomes and residuals use the observed structural design `X`. This distinction
matters: ordinary OLS software applied to first-stage fitted exposures can
reproduce the 2SLS coefficient, but its default residual variance and standard
errors are not valid IV inference.

## Check instrument strength and identification

Each endogenous column has an excluded-instrument diagnostic:

```java
var strength = fit.strengthDiagnostics().get(0);
System.out.printf(
    "partial R2=%.4f, classical F=%.2f, covariance-matched F=%.2f, p=%.4g%n",
    strength.partialRSquared(),
    strength.classicalFStatistic(),
    strength.effectiveFStatistic(),
    strength.pValue());

System.out.println("minimum first-stage F = "
    + fit.minimumFirstStageFStatistic());
System.out.println("projected-design condition number = "
    + fit.projectedDesignConditionNumber());
```

The partial R-squared compares the full first stage to the included exogenous
regressors alone. `classicalFStatistic()` is the conventional homoskedastic
excluded-instrument F. `effectiveFStatistic()` instead matches the selected
covariance: it is the HC Wald chi-square divided by its numerator degrees of
freedom for `HC0`/`HC1`, and a finite-cluster Wald F for clustered fits. Inspect
`referenceDistribution()` and `denominatorDegreesOfFreedom()` before reporting
the test.

The familiar first-stage F greater than 10 rule is only a rough screen. It is
not a universal weak-instrument guarantee, particularly with multiple
endogenous regressors, heteroskedasticity, or clustering. This API reports
first-stage diagnostics and numerical rank, but it does not provide
weak-identification-robust Anderson-Rubin/CLR confidence sets. If weak
identification is plausible, use a method designed for it rather than treating
a conventional 2SLS Wald interval as conclusive.

JLinAlg rejects all of these cases before returning a fit:

- fewer excluded instruments than endogenous regressors;
- a rank-deficient structural or instrument matrix;
- a rank-deficient projected structural design, meaning the instruments do not
  identify every structural coefficient; or
- too few observations for structural or first-stage residual degrees of
  freedom.

`projectedDesignConditionNumber()` scales each projected column by the norm of
its observed structural counterpart, so it is invariant to simple unit changes
while retaining information about a weak projection. A very large finite value
warns that identification is numerically fragile even when the rank check
passes.

## Cluster-robust inference

Supply one integer cluster identifier per row when observations may be
dependent within groups:

```java
int[] household = {0, 0, 1, 1, 2, 3, 4, 5};
InstrumentalVariableOptions clusteredOptions =
    new InstrumentalVariableOptions(
        InstrumentalVariableCovariance.CLUSTER_CR1, 0.95);

InstrumentalVariableResult clustered = InstrumentalVariableRegression.fit(
    outcome, exogenous, endogenous, excludedInstruments,
    household, clusteredOptions, BackendPolicy.CPU);
```

`CLUSTER_CR0` accumulates the outer product of cluster-level IV scores.
`CLUSTER_CR1` additionally multiplies it by
`G/(G-1) * (n-1)/(n-p)`. Coefficient tests use `G-p` denominator degrees of
freedom, matching the library's cluster-t convention. First-stage excluded-
instrument tests use `G-l`, where `l` is the number of included plus excluded
instrument columns. Consequently, clustered fitting requires more independent
clusters than both the structural and first-stage parameter counts. Clusters
must be independent; the correction does not rescue a design with only a few
dominant clusters.

## Multiple endogenous regressors and contrasts

Put every endogenous regressor in a separate column. There must be at least as
many excluded instrument columns as endogenous columns. The result exposes the
first-stage fitted matrix and a diagnostic for each endogenous column:

```java
double[][] firstStageFits = fit.instrumentedEndogenous();

// Jointly test the exogenous slope and endogenous exposure coefficient.
var joint = fit.testContrast(new double[][] {
    {0, 1, 0},
    {0, 0, 1}
});
```

Coefficient-level tests use asymptotic normal inference for HC0/HC1,
residual-t inference for homoskedastic fits, and cluster-t inference for
clustered fits. Joint `testContrast` calls use a Wald chi-square law for HC0/HC1
and an F law for homoskedastic or clustered fits. First-stage column-wise F
statistics do not replace a substantive assessment of joint identification
when several endogenous regressors are present.

## Reproduce the independent R fixture

The checked-in fixture uses two endogenous regressors, three excluded
instruments, 30 clusters, and deterministic heteroskedastic errors. Base R 4.6.1
computes the pivoted-QR projection and explicit IV/GMM sandwiches independently
for all five covariance choices:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/test/resources/r-reference/generate-instrumental-variable-reference.R

.\gradlew.bat test `
  --tests org.jlinalg.regression.InstrumentalVariableRegressionTest `
  --tests org.jlinalg.regression.InstrumentalVariableRReferenceTest
```

The generator requires only base R and `stats`. Its outputs record the R runtime
and calculation contract alongside the numeric fixture.

## Scope

This API is linear 2SLS with finite rectangular, complete data. It does not yet
implement formula parsing, observation weights, missing-data omission, LIML,
continuously updated GMM, nonlinear endogenous models, weak-IV-robust tests, or
instrument exogeneity tests such as Hansen's J. Scale and encode inputs before
fitting, and report the covariance choice, cluster count where applicable,
first-stage diagnostics, and identifying assumptions with the estimate.
