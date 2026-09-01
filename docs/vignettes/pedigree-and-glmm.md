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

The variance order is pedigree terms, ordinary terms, then residual. The
sparse combined result currently does not provide the dense animal model's
full scalable PEV/reliability matrix.

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

Satterthwaite and Kenward-Roger options describe the final PQL working model,
not an exact finite-sample distribution for the original non-Gaussian model.
