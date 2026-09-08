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

Stopping is checked in fitted-response space, rather than raw parameter units:
the projected Gauss–Newton correction must satisfy both the parameter tolerance
and the quadratic model's relative-objective tolerance, relative to the residual
norm. There is no absolute response-magnitude floor: adding a large constant to
the response and mean must not certify an unfitted model. Fixed-model
backtracking requires a strict SSE decrease. A rounding-limited plateau retains
the best accepted parameters but reports `converged() == false` unless the
stationarity criterion is satisfied. Center known offsets before fitting when
possible; at offsets such as `1e15`, evaluating the mean itself can quantize
meaningful residual changes. Tiny changes after backtracking alone do not
certify a solution.
Inference uses the Jacobian evaluated at the **returned** parameters, including
when the iteration budget expires. Fixed fits require more observations than
parameters so residual-variance estimation has positive degrees of freedom.

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

Backtracking honors `initialStep` and `maximumStepHalvings`, comparing the
nonlinear conditional mean at a fixed current BLUP. Raw residual sums from
different variance/BLUP fits are not used as a common likelihood objective.
The final linearization is refreshed after the final accepted parameter change;
`converged()` requires both the nonlinear stationarity check and the underlying
variance fit to converge. `objective()` is the conditional residual sum of
squares, not a marginal log likelihood. Sparse random-effect contributions are
recovered from the linearized fitted vector without expanding `Z` to a dense
observation-by-random-coefficient matrix.

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

## Independent audit evidence

The September 2026 audit includes noisy-exponential `stats::nls` coefficient,
standard-error and SSE checks; an independent dense Gaussian ML calculation
for a nonlinear mean with correlated pedigree effects; scaling and iteration
budget counterexamples; and mixed-step-control regression tests. Pedigree ML
parameters differ from the independent R optimizer by at most `2.3e-8`, and
the fitted log likelihood agrees to about `4e-14` on that fixture. This is
fixture-specific validation of the additive model, not general `nlme` parity.

Reproduction commands, warm timings, checksums, and limitations are in the
[audit evidence](../../src/benchmark/resources/remaining-model-audit-benchmark/audit-results.md).
