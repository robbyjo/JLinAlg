# Pedigree models and generalized mixed models

## Build and validate a pedigree

Individuals may arrive in any order. Parents must either be present or `null`:

```java
Pedigree pedigree = Pedigree.of(List.of(
    PedigreeIndividual.founder("sire"),
    PedigreeIndividual.founder("dam"),
    new PedigreeIndividual("offspring", "sire", "dam")));

double relationship = pedigree.relationship("sire", "offspring");
double[] inbreeding = pedigree.inbreedingCoefficients();
SparseSymmetricMatrix aInverse =
    pedigree.sparseRelationshipMatrixInverse();
```

The constructor rejects duplicate identifiers, unknown named parents, and
ancestry cycles. Unknown founder sources are represented with `null`, not a
string sentinel such as `"0"`.

## Dense animal-model REML

Observation identifiers define the incidence matrix and may repeat. Pedigree
members without phenotypes remain eligible for breeding-value prediction:

```java
double[] y = {8.5, 9.1, 8.9, 9.5, 10.8, 11.4};
double[][] x = {{1}, {1}, {1}, {1}, {1}, {1}};
List<String> animal = List.of(
    "sire", "sire", "dam", "dam", "offspring", "offspring");

PedigreeRemlResult fit = PedigreeReml.fit(y, x, animal, pedigree);
double h2 = fit.heritability();
double offspringBlup = fit.breedingValue("offspring");
double offspringReliability = fit.reliability("offspring");
Map<String, Double> animalModes = fit.ranef();
```

The dense reference path returns additive/residual variance, heritability,
BLUP, PEV, and reliability. Check the underlying REML convergence before using
the estimates scientifically.

Retain the animal-model structure for simulation and bootstrap:

```java
PreparedPedigreeReml prepared = new PreparedPedigreeReml(
    x, animal, pedigree, RemlOptions.defaults(),
    BackendPolicy.PREFERRED);
PedigreeRemlResult fitted = prepared.fit(y);

double[][] simulated = PedigreeSimulation.simulate(
    prepared, fitted, 100,
    42, MixedModelSimulationMode.MARGINAL);
GaussianBootstrapResult intervals = PedigreeBootstrap.bootstrap(
    prepared, fitted,
    new BootstrapOptions(999, 0.95, 42, 1));
```

Marginal simulation draws breeding values jointly from the fitted numerator-
relationship covariance. Conditional simulation retains the BLUPs. Pedigree
results also expose `fixef()`, named `ranef()`, `varCorr()`, and conditional
`fittedValues()`/`residuals()` aliases.

## Sparse animal model and additional random terms

For large pedigrees, variance estimation can consume `A^-1` directly:

```java
SparsePedigreeRemlResult sparse = SparsePedigreeReml.fit(
    y, x, animal, pedigree,
    RemlOptions.builder().initialVariances(1, 1).build(),
    BackendPolicy.PREFERRED);

double[] pev = sparse.predictionErrorVariances();
double[] reliability = sparse.reliabilities();
Map<String, Double> selected = sparse.predictionErrorVariances(
    List.of("member-17", "member-42"));
```

To combine several pedigree structures with ordinary effects, construct named
pedigree terms and fit one sparse precision model:

```java
PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
    "animal", animal, pedigree);
RandomEffectTerm batch = RandomEffectTerm.randomIntercept(
    "batch", List.of("a", "a", "b", "b", "c", "c"));

SparseLinearMixedModelResult combined = SparsePedigreeMixedModel.fit(
    y, x, List.of(additive), List.of(batch),
    RemlOptions.builder().initialVariances(1, 1, 1).build(),
    BackendPolicy.PREFERRED);
```

The variance order is pedigree terms, ordinary terms, then residual. Sparse
diagonal PEVs are extracted in fixed-size solve batches without materializing
the full inverse, so working memory is linear in coefficient count times the
batch size. A complete covariance matrix remains intentionally unavailable
because its storage is quadratic; request selected individuals when only
named diagonal uncertainty is needed.

## Generalized linear mixed models

`GlmmPql` handles binomial, Poisson, negative-binomial, and supported
quasi-family mixed models by first-order PQL:

```java
double[] counts = {1, 2, 1, 5, 4, 6};
double[][] fixed = {{1}, {1}, {1}, {1}, {1}, {1}};
double[] groupRelationship = {
    1, 1, 1, 0, 0, 0,
    1, 1, 1, 0, 0, 0,
    1, 1, 1, 0, 0, 0,
    0, 0, 0, 1, 1, 1,
    0, 0, 0, 1, 1, 1,
    0, 0, 0, 1, 1, 1
};

GlmmPqlResult pql = GlmmPql.fit(
    counts, fixed, GlmFamilies.poisson(),
    List.of(new VarianceComponent(
        "group", counts.length, groupRelationship)));
```

Offsets and prior weights are available in the full overload. PQL inference is
conditional on the final working Gaussian linearization. It is not a Laplace
or adaptive-quadrature marginal GLMM likelihood and can be biased for rare
binary outcomes or small clusters.

Use pedigree covariance in the same PQL engine with:

```java
PedigreeGlmmPqlResult pedigreePql = PedigreeGlmmPql.fit(
    binaryResponse, fixed, GlmFamilies.binomial(),
    animal, pedigree);
```

For a first-order Laplace marginal likelihood, build grouped designs and
coefficient-space precision matrices directly:

```java
RandomEffectTerm levySet = RandomEffectTerm.randomIntercept(
    "Levy_Set", levySetIds);
PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
    "animal", animal, pedigree);

try (SparseGlmmLaplace.Prepared scan =
        SparseGlmmLaplace.prepareWithPrecision(
            binaryResponse.length,
            GlmFamilies.binomial(),
            List.of(levySet, additive.randomEffect()),
            List.of(SparsePrecisionMatrix.identity(
                levySet.coefficients()), additive.precision()),
            GlmmLaplaceOptions.defaults(),
            BackendPolicy.PREFERRED)) {
    GlmmLaplaceResult fit = scan.fit(
        binaryResponse, rowMajorFixedDesign, fixedColumns);
}
```

`Prepared` owns one backend for the scan and lazily creates one reusable sparse
Cholesky factor per calling worker. The random-effect Hessian is assembled as
`Z'WZ + Q`, and fixed effects are solved through its Schur complement. Pedigree
`A^-1` and grouped `Z` therefore remain sparse; neither `A` nor `ZAZ'` is
materialized. Close the prepared scan only after all worker tasks finish.

## Adaptive quadrature for independent grouped random intercepts

`GlmmQuadrature` maximizes the marginal likelihood with mode/curvature adaptive
Gauss–Hermite integration. It accepts the built-in `GlmFamilies.binomial()`
(logit, Bernoulli or binomial counts) and `GlmFamilies.poisson()` (log, counts),
both at dispersion one. Other families, including quasi-likelihood and families
with unknown dispersion, are rejected. Groups must be independent; fitting
requires at least two groups and a full-rank fixed design.

```java
GlmmQuadratureOptions controls = GlmmQuadratureOptions.defaults();
GlmmQuadratureResult fit = GlmmQuadrature.fit(
    binaryResponse, fixed, groupIds, GlmFamilies.binomial(), controls);
if (!fit.converged()) throw new IllegalStateException(fit.status());
double[] beta = fit.beta();
double[] se = fit.standardErrors();
double[][] covariance = fit.fixedEffectCovariance();
double variance = fit.randomVariance();
```

The full overload accepts trial counts and offsets before `controls`. Binomial
responses are proportions with integer `response * trials` (allowing two ULPs
of round-trip floating-point error for supplied trials greater than one); without trials,
responses must be zero or one. A binomial count row includes its combinatorial
constant, so collapsing ordered Bernoulli rows changes the reported likelihood
by the corresponding log binomial coefficient. Poisson trials must be absent
or all one; exposure belongs in the log offset. General observation weights,
random slopes, crossed grouping factors, and estimated dispersion are not
implemented in this API.

For each group, write `u = b / sd` and
`h(u) = sum(log p(y | X beta + offset + sd*u)) - u*u/2`.
Safeguarded Newton iterations locate its unique mode `m`; the scale is
`a = 1 / sqrt(-h''(m))`. The marginal integral is approximated by
`a / sqrt(pi) * sum(w[k] * exp(h(m + sqrt(2)*a*z[k]) + z[k]*z[k]))`.
Newton proposals must contract the score bracket, so large Poisson offsets do
not consume the iteration budget in unit-sized steps. The implementation uses
stable log probabilities, compensated sums, and log-sum-exp. Large-count masses
use Stirling-error normalizers and cancellation-safe deviance terms instead of
subtracting large log factorials. Hermite rules come from a tridiagonal eigensolver
and are cached.

Defaults start at 9 nodes and refine to at most 257, requiring two successive
agreements within a per-group allocation of `1e-10` absolute total
log-likelihood tolerance. `GlmmQuadratureOptions(initialNodes, maximumNodes,
quadratureTolerance, maximumIterations, gradientTolerance)` controls these
limits; the maximum supported order is 512. `evaluate(...)` exposes the
fixed-parameter likelihood, maximum order actually used, estimated error, and
convergence. A one-node rule reproduces Laplace, but cannot establish integration
accuracy. Node-refinement error is an estimate, not a rigorous bound. Exhausting
the budget never counts as quadrature convergence. The zero-variance model is
evaluated directly without quadrature.

Fitting uses projected BFGS with design- and information-scaled coefficients,
two interior variance starts, and an explicitly optimized zero-variance candidate.
The variance coordinate is `t = log1p(variance / v0)`, where `v0` is the reciprocal
of the largest group's reference information. This resolves concentrated
optima without imposing a positive variance floor. Interior derivatives
differentiate the full marginal likelihood numerically. At variance zero,
the coefficient score and right variance derivative are analytic; the latter
uses `d log L / d variance = sum(h'' + h'^2) / 2`, with derivatives of each
group's conditional log likelihood with respect to its random intercept.
Convergence requires a small
projected score, small information-based coefficient correction, successful
integration, and positive definite observed information. Exhausted iterations,
line searches, separation, and singular information do not automatically become
successful fits. Multiple starts improve robustness but do not certify a global
maximum. The iteration limit applies to each optimization start.

At an interior optimum, `parameterCovariance()` transforms the inverse numerical
Hessian back to the order fixed coefficients, random **variance**. Its fixed-effect block
accounts for nuisance-variance uncertainty. `waldZ()` and `pValues()` provide
asymptotic normal fixed-effect tests. At exactly zero variance,
`varianceBoundary()` is true and `jointInferenceAvailable()` is false: only the
fixed-effect covariance conditional on variance zero is returned, and the
variance row/column is `NaN`. No ordinary Wald variance test is supplied at that
nonregular boundary. Nonconverged estimates have `NaN` covariance.

`fittedMeans()` averages over a **new** group's random-intercept distribution;
it returns binomial probabilities or Poisson means, rather than conditional
posterior-mode predictions. Binomial prediction integrals are separately checked
and return `NaN` if they exhaust their node budget. Poisson means use the exact
lognormal expectation `exp(eta + variance/2)`.

The checked-in R fixtures cover eight fixed-parameter integrals, including
10,000 Bernoulli observations, rare/all-zero outcomes, and concentrated Poisson
counts, plus four 384-observation fits against `lme4::glmer(nAGQ=25)`. They check
fixed effects, variance, full likelihood, joint covariance, standard errors,
and new-group means. The audit's 50/50 Bernoulli example now gives
`-70.94145611404414`, agreeing with R numerical integration. For count data,
`lme4`'s reported nAGQ>1 likelihood uses a different additive constant; references
independently integrate the full binomial/Poisson probability mass.

Regenerate with `src/test/resources/r-reference/generate-quadrature-reference.R`.
The additional `generate-quadrature-stress-reference.R` uses base R to cover
both review-reported concentrated interior optima, their likelihood and covariance,
count masses through `1e16`, valid large-denominator proportions, and Poisson
offsets through 1000. Regression tests also reject genuinely fractional counts.
Run `org.jlinalg.benchmark.GlmmQuadratureBenchmark` and
`src/benchmark/r/quadrature_benchmark.R` from the repository root for comparable
full-fit timings with accuracy gates. The recorded run and limitations are in
`src/benchmark/resources/quadrature/accuracy-and-timing.md`.

This is a numerical adaptive-quadrature alternative for independent scalar
random intercepts. It does **not** implement pedigree-correlated or
multidimensional adaptive quadrature. Pedigree models in this vignette still use
the separate PQL or sparse Laplace paths described above; their approximation
limitations are not removed by this addition.

## Zero-inflated Poisson and negative-binomial mixed models

`SparseZeroInflatedMixedModel` maximizes a frequentist first-order Laplace
likelihood. BOBYQA sees only fixed effects, NB2 size coefficients, and
log-variance components. Each objective evaluation finds the potentially large
random-effect mode with sparse damped Newton iterations; its observed Hessian
is reused for the Laplace determinant rather than computing expected Fisher
information for BOBYQA.

Independent ordinary or pedigree random effects may enter either predictor.
For a correlated two-process NB2 pedigree model:

```java
PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
    "individual", observationIndividualIds, pedigree);

CorrelatedZeroInflatedRandomEffect joint =
    CorrelatedZeroInflatedRandomEffect.pedigree("additive", additive);

try (SparseZeroInflatedMixedModel.Prepared prepared =
        SparseZeroInflatedMixedModel.prepareNegativeBinomial(
            count.length, List.of(), List.of(), List.of(), List.of(),
            List.of(joint), ZeroInflatedMixedOptions.defaults(),
            BackendPolicy.PREFERRED)) {
    ZeroInflatedMixedResult fit = prepared.fitWithInference(
        count, countFixed, countColumns,
        zeroFixed, zeroColumns, sizeFixed, sizeColumns, countOffset);
}
```

The NB2 conditional variance is `mu + mu^2 / size`. The result distinguishes
the conditional count mean, structural-zero probability, unconditional fitted
mean, and total fitted zero mass. Pedigree `A^-1` remains in coefficient space,
so unobserved ancestors are retained without forming dense `A` or `ZAZ'`.

For correlated effects, `G` has separate count/zero variances and a correlation
represented as `0.99 * tanh(z)`. The random precision is
`inverse(G) kron inverse(A)` and stays sparse. Zero observations contribute the
full observed 2-by-2 predictor curvature, including the cross term; positive
observations have a zero cross term.

Automatic outer optimization uses bounded BFGS and workload-gated parallel
numerical gradients for modes with at most 128 random coefficients. Each
gradient perturbs the complete Laplace objective from the same converged center
mode, so it includes mode and log-determinant changes. Larger sparse modes use
BOBYQA; callers may also select BOBYQA explicitly or limit gradient threads in
`ZeroInflatedMixedOptions`. `fitWithInference` performs a separate post-fit
numerical Hessian, including the variance/correlation nuisance parameters, and
exposes fixed/dispersion covariance and standard errors. `profile` refits the
nuisance parameters over a caller grid; `parametricBootstrap` records failed
replicates as `NaN`. Prepared state holds one symbolic analysis and lazily
creates one numeric factor per calling worker.

`ZeroInflatedMixedFormula` compiles count and zero formulas, plus an NB2 size
formula when needed. Checked-in fixtures compare grouped ZIP/ZINB fits and ZIP
standard errors with `glmmTMB`; a separately compiled TMB sparse-GMRF template
checks correlated pedigree fixed effects, variances, correlation, and marginal
likelihood. The default benchmark fits 1,500 observations and 1,000 random
coefficients with 4,900 sparse equation/factor entries in 2.161532 seconds on
the documented development host; use `benchmarkZeroInflatedMixed` to rerun it.

Satterthwaite and Kenward-Roger options describe the final PQL working model,
not an exact finite-sample distribution for the original non-Gaussian model.
