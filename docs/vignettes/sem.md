# Structural equation modeling

JLinAlg 0.2.0 fits continuous observed-variable path models by Gaussian
covariance-structure maximum likelihood. The engine uses the RAM relationship

`Sigma = (I - A)^-1 S (I - A)^-T`,

where `A` contains directed paths and `S` contains residual variances and
covariances. This vignette builds a model, fits rows or a covariance matrix,
checks inference and fit, and explains the current scope.

## Specify a path model

Variable order defines both input columns and covariance-matrix order:

```java
SemModel model = SemModel.builder(
        "Age", "Sex", "BMI", "Waist", "Glucose")
    .regression("BMI", "Age", 0.1)
    .regression("BMI", "Sex", -0.1)
    .regression("Waist", "BMI", 0.7)
    .regression("Waist", "Sex", -0.1)
    .regression("Glucose", "Age", 0.1)
    .regression("Glucose", "BMI", 0.2)
    .covariance("Age", "Sex", 0.0)
    .variance("Age", 1.0)
    .variance("Sex", 1.0)
    .variance("BMI", 0.8)
    .variance("Waist", 0.3)
    .variance("Glucose", 0.8)
    .build();
```

The numeric values are optimizer starts. Use `fixedRegression`,
`fixedVariance`, or `fixedCovariance` when a value is known rather than
estimated. When a variance is omitted, the builder adds a free variance with
start 1.0.

Repeated labels impose equality constraints:

```java
SemModel equalPaths = SemModel.builder("x", "m", "y")
    .regression("shared", "m", "x", 0.4)
    .regression("shared", "y", "m", 0.4)
    .build();
```

## Fit observed rows

Rows are observations and columns follow the model variable order:

```java
SemFitResult fit = Sem.fit(
    observedRows,
    model,
    new SemOptions(10_000, 1e-8, MissingDataPolicy.OMIT),
    BackendPolicy.PREFERRED);

if (!fit.converged()) {
    throw new IllegalStateException("SEM did not converge");
}
```

`MissingDataPolicy.OMIT` uses one complete-case sample across every model
variable. It is not full-information maximum likelihood. Standardize variables
when scales differ greatly; this improves optimizer conditioning without
changing standardized substantive relationships.

The row-data likelihood centers each selected column and forms its ML
covariance with denominator `N`. Means are not modeled.

## Fit a covariance matrix

Use `fitCovariance` when sufficient statistics are already available:

```java
double[] covariance = {
    1.00, 0.30, 0.25,
    0.30, 1.00, 0.55,
    0.25, 0.55, 1.00
};

SemModel mediation = SemModel.builder("x", "m", "y")
    .regression("m", "x", 0.3)
    .regression("y", "m", 0.5)
    .build();

SemFitResult covarianceFit = Sem.fitCovariance(
    covariance, 2_000, mediation,
    SemOptions.defaults(), BackendPolicy.PREFERRED);
```

The matrix is dense row-major and must be positive definite. Pass the ML
covariance, not an unadjusted `N - 1` covariance, when exact row-data
likelihood and fit-index parity is required.

## Inspect parameters and fit

```java
SemParameterEstimate bmiToWaist = fit.parameter("Waist~BMI");
System.out.printf(
    "beta=%g se=%g z=%g p=%g%n",
    bmiToWaist.estimate(),
    bmiToWaist.standardError(),
    bmiToWaist.zStatistic(),
    bmiToWaist.pValue());

System.out.printf(
    "chi-square=%g df=%d CFI=%g TLI=%g RMSEA=%g SRMR=%g%n",
    fit.chiSquare(), fit.degreesOfFreedom(),
    fit.cfi(), fit.tli(), fit.rmsea(), fit.srmr());
```

Results also expose the implied covariance, log likelihood, AIC, BIC, sample
size, function-evaluation count, and selected compute backend. A
just-identified model has zero degrees of freedom and can reproduce the sample
covariance by construction; its global fit indices do not test a restrictive
model.

## TOPMed worked benchmark

The reproducible benchmark uses 4,680 complete observations on `Sex`, `Age`,
`BMI`, `Waist`, `Systolic_BP`, `Diastolic_BP`, `Glucose`, `HDL`, `LnTG`,
`LogInsulin`, `CRP`, and `eGFR`. The full requested 57-variable set has only
five jointly complete rows and is not estimable by the current complete-case
engine.

For the 36-parameter model, JLinAlg and `lavaan` 0.7-2 produced identical
printed log likelihood (-70763.3172885), chi-square (4497.68284338), CFI
(0.799908468824), TLI (0.685570451009), and RMSEA (0.150559985169).
Maximum absolute differences were below `1e-8` for estimates and `5e-10` for
standard errors. Median warmed fit times were 0.0081623 seconds for JLinAlg
and 0.0700000 seconds for `lavaan`, an 8.58x speedup on the documented host.

The poor CFI and RMSEA are a substantive warning about this illustrative path
structure, not a numerical discrepancy: both engines found the same optimum.
See [TOPMed SEM validation and performance](../topmed-sem-performance.md) for
the complete model, results, environment, and reproduction commands.

## Numerical validation

`SemTest` includes a deterministic `lavaan` 0.7-2 fixture with directed paths,
free variances, an exogenous covariance, a correlated disturbance,
expected-information standard errors, likelihood, and fit indices. This
guards the general RAM derivative and inference path independently of the
TOPMed example.

## Scope and limitations

The current engine supports continuous observed variables, directed paths,
free or fixed variances and covariances, equality labels, complete-case
covariance ML, Wald inference, and conventional global fit indices.

It does not implement latent measurement variables, intercept/mean structures,
ordinal thresholds, robust or clustered corrections, modification indices,
indirect-effect delta-method inference, or FIML missing-data patterns. Binary
variables are treated as Gaussian continuous observations. Use a package with
the required likelihood rather than interpreting these unsupported cases as
silently approximated.
