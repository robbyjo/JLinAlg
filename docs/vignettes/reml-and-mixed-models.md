# REML and Gaussian mixed models

## General covariance-component REML

REML models `V` as a nonnegative combination of known covariance bases. This
example uses a grouping relationship plus independent residual variance:

```java
double[] y = {0, 1, 2, 4, 5, 6, 8, 9, 10};
double[][] x = {
    {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
};
List<String> subject = List.of(
    "a", "a", "a", "b", "b", "b", "c", "c", "c");

List<VarianceComponent> components = List.of(
    VarianceComponent.randomIntercept("subject", subject),
    VarianceComponent.identity("residual", y.length));

RemlResult fit = Reml.fit(y, x, components);
if (!fit.converged()) {
    throw new IllegalStateException(fit.convergenceMessage());
}
System.out.println(Arrays.toString(fit.varianceComponents()));
System.out.println(Arrays.toString(fit.beta()));
```

Caller-supplied positive-semidefinite bases support kinship, repeated-measure,
spatial, or other scientifically specified covariance structures. Do not add
an identity basis twice.

## Degrees of freedom and ML

The fast default is `N-rank(X)-1`. Small samples can request Satterthwaite or
Kenward-Roger inference:

```java
RemlOptions smallSample = RemlOptions.builder()
    .degreesOfFreedomMethod(DegreesOfFreedomMethod.KENWARD_ROGER)
    .build();
RemlResult adjusted = Reml.fit(
    y, x, components, smallSample, BackendPolicy.PREFERRED);
```

Use `DegreesOfFreedomMethod.SATTERTHWAITE` for its delta-method alternative.
Both are slower and dense. For likelihood-ratio comparison of nested fixed
effects, fit both models with `VarianceEstimation.ML`; REML likelihoods from
different fixed-effect designs are not comparable.

## Term-oriented LMM

`LinearMixedModel` constructs covariance bases and conditional modes from
random-effect terms:

```java
double[] age = {0, 1, 2, 0, 1, 2, 0, 1, 2};
List<RandomEffectTerm> random = List.of(
    RandomEffectTerm.randomIntercept("subject", subject),
    RandomEffectTerm.randomSlope("subject:age", subject, age));

LinearMixedModelResult lmm = LinearMixedModel.fit(y, x, random);
double[] subjectModes = lmm.randomEffects("subject").estimates();
double[] conditionalResiduals = lmm.residuals();
```

These two terms have independent variances. When an intercept/slope covariance
is required, construct a `CorrelatedRandomEffectBlock` and fit it with
`CorrelatedLinearMixedModel`. The correlated-block reference likelihood is
currently dense.

For many groups and independent terms, use `SparseLinearMixedModel`:

```java
SparseLinearMixedModelResult sparse = SparseLinearMixedModel.fit(
    y, x, random,
    RemlOptions.builder().initialVariances(1, 1, 1).build(),
    BackendPolicy.PREFERRED);
System.out.println(sparse.equationNonzeroCount());
```

The initial-variance count is one per random term plus the residual. The sparse
fitter currently uses the fast residual-DF approximation; use the dense path
when Satterthwaite, Kenward-Roger, or full PEVs are required.

## Prediction on new data

Supply the fixed design and random terms in the same coefficient naming
convention as the fit:

```java
double[][] newX = {{1}, {1}};
RandomEffectTerm newSubjects = RandomEffectTerm.randomIntercept(
    "subject", List.of("a", "previously-unseen"));

double[] population = MixedModelPrediction.marginal(lmm, newX);
double[] conditional = MixedModelPrediction.conditional(
    lmm, newX, List.of(newSubjects), true);
```

With `allowNewLevels=true`, an unseen level gets random mode zero, so its
conditional prediction equals its population prediction. Set it to `false` to
fail on unseen levels. Omitting a fitted random term from the prediction list
is equivalent to excluding that term from the conditional prediction.

## Refit many responses

Retain the design once and warm-start new responses from a previous fit:

```java
PreparedLinearMixedModel prepared = new PreparedLinearMixedModel(
    x, random, RemlOptions.defaults(), BackendPolicy.PREFERRED);

LinearMixedModelResult first = prepared.fit(y);
LinearMixedModelResult next = prepared.refit(first, anotherResponse);
```

This is useful for repeated phenotypes with identical rows and model structure.
For large association scans, the dedicated association engines reuse still
more null-model work; see the [GWAS/TWAS vignette](association-gwas-twas.md).

## Simulate and bootstrap a fitted model

Marginal simulation draws new random effects and residuals. Conditional
simulation retains fitted random-effect modes and draws only residuals:

```java
double[][] marginal = MixedModelSimulation.simulate(
    prepared, first, 100, 42,
    MixedModelSimulationMode.MARGINAL);
double[][] conditional = MixedModelSimulation.simulate(
    prepared, first, 100, 42,
    MixedModelSimulationMode.CONDITIONAL);
```

Parametric bootstrap simulates marginally, warm-refits the prepared model, and
summarizes successful replicates:

```java
BootstrapOptions bootstrap = new BootstrapOptions(
    999, 0.95, 20260901L, 8);
GaussianBootstrapResult intervals = MixedModelBootstrap.bootstrap(
    prepared, first, bootstrap);

BootstrapParameterSummary slope =
    intervals.fixedEffectSummaries().get(1);
```

Each summary reports the fitted estimate, bootstrap mean, bias, empirical SE,
and percentile interval. Failed/nonconverged replicates remain in
`failures()`. Use one bootstrap worker with GPU/native multithreaded BLAS; raise
parallelism deliberately for a CPU backend.
