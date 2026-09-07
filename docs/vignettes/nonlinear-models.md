# Nonlinear fixed and mixed effects

JLinAlg provides deterministic Gaussian nonlinear least squares through an
analytic-gradient Gauss–Newton solver. The same mean function can be fitted
with sparse ordinary random effects or sparse pedigree effects. The mixed
solver linearizes the nonlinear fixed mean at each iteration and reuses a
prepared sparse REML equation structure for the random effects.

## Nonlinear fixed effects

`NonlinearMeanFunction` receives the current parameter vector and row index and
returns both the mean and its analytic gradient:

```java
NonlinearMeanFunction exponential = (theta, row) -> {
    double x = row * 0.25;
    double mean = theta[0] * Math.exp(theta[1] * x);
    return new NonlinearMeanFunction.Evaluation(mean,
        new double[] {Math.exp(theta[1] * x), x * mean});
};

NonlinearFitResult fit = NonlinearFixedModel.fit(
    response, new double[] {1.0, 0.1}, exponential,
    NonlinearModelOptions.defaults(), BackendPolicy.PREFERRED);
```

The result exposes nonlinear parameter estimates through `beta()`, local
linearized standard errors, fitted values, residuals, objective, and
convergence metadata. Analytic gradients are required so rank and convergence
failures are visible instead of being hidden behind finite-difference noise.

## Ordinary nonlinear mixed effects

Add a sparse random-effect term to the same mean function:

```java
RandomEffectTerm site = RandomEffectTerm.randomIntercept("1|site", siteIds);
NonlinearMixedFitResult mixed = NonlinearMixedModel.fit(
    response, new double[] {1.0, 0.1}, exponential,
    List.of(site), RemlOptions.defaults(),
    NonlinearModelOptions.defaults(), BackendPolicy.PREFERRED);
```

The random effects remain additive on the response scale:

```text
y = f(theta, x) + Z b + e
```

The `linearizedModel()` result retains sparse variance components, conditional
modes, and sparse diagnostics from the final linearization.

## Pedigree nonlinear mixed effects

Pedigree effects use the existing sparse numerator-relationship precision and
can be combined with ordinary random terms:

```java
PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
    "additive", observationIds, pedigree);
NonlinearMixedFitResult animal = NonlinearMixedModel.fitPedigree(
    response, new double[] {1.0, 0.1}, exponential,
    List.of(additive), List.of(), RemlOptions.defaults(),
    NonlinearModelOptions.defaults(), BackendPolicy.PREFERRED);
```

Unphenotyped ancestors remain in the sparse coefficient system, so their
relationship information contributes to the observed-animal fit without
materializing a dense observation covariance matrix.

## Statistical boundary

This is a Gaussian additive nonlinear mixed-effects path using first-order
Gauss–Newton linearization and REML variance updates. It is not a Laplace or
adaptive-Gauss–Hermite nonlinear likelihood for non-Gaussian outcomes, and it
does not claim exact `nlme` or `lme4` nonlinear likelihood parity. Use the
convergence flag, inspect the final linearized model, and validate the
parameterization against an independent reference fixture for each model.
