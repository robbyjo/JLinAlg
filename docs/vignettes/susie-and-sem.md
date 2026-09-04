# SuSiE fine mapping and structural equation models

> **v0.1.0 performance status:** These implementations are tested, but neither
> SuSiE nor SEM is yet optimized or benchmarked for large or repeated
> workloads.

## SuSiE with individual-level data

Rows are samples and columns are candidate variables in the same order as the
names list:

```java
SusieOptions options = new SusieOptions(
    5,      // maximum single effects
    200,    // iterations
    1e-8,   // convergence tolerance
    0.2,    // prior effect variance
    true,   // estimate residual variance
    0.95,   // credible-set coverage
    0.5);   // minimum credible-set purity

SusieResult fineMap = Susie.fit(
    phenotype, genotypeMatrix, variantNames,
    options, BackendPolicy.PREFERRED);

if (!fineMap.converged()) {
    throw new IllegalStateException("SuSiE did not converge");
}
double[] pip = fineMap.pip();
double[] posteriorBeta = fineMap.posteriorMean();
List<CredibleSet> sets = fineMap.credibleSets();
```

Individual-level inputs are centered and scaled once. Posterior means are
returned on the original predictor scale with a fitted intercept.

## Summary z scores and LD

```java
double[] z = {12.0, 0.2, -0.1, -10.0, 0.3};
double[][] ld = {
    {1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0},
    {0, 0, 1, 0, 0},
    {0, 0, 0, 1, 0},
    {0, 0, 0, 0, 1}
};

SusieResult summary = Susie.fitSummary(
    z, ld, 1_000, List.of("a", "b", "c", "d", "e"),
    new SusieOptions(2, 200, 1e-8, 0.2, false, 0.95, 0.5),
    BackendPolicy.PREFERRED);
```

Summary mode assumes standardized quantitative traits/predictors and a
positive-semidefinite LD correlation matrix from an ancestry-matched reference.
Use `fitSufficientStatistics` when `X'X`, `X'y`, and `y'y` are available on a
different scale. Credible-set purity is the minimum absolute within-set LD;
always inspect it alongside coverage and PIP.

## Colocalize SuSiE signals

Two JLinAlg fits can be colocalized directly:

```java
ColocSusieResult coloc = ColocSusie.analyze(exposureFineMap, outcomeFineMap);

for (ColocSignalPair pair : coloc.signalPairs()) {
    System.out.printf("effects %d/%d: H4=%.4f%n",
        pair.trait1EffectIndex(), pair.trait2EffectIndex(),
        pair.posteriorH4());
}
double[] sharedVariantPosterior = coloc.sharedVariantPosterior(0);
```

`ColocSusieInput` also accepts an L-by-P matrix of log Bayes factors, which is
the direct equivalent of selected rows from `susieR`'s `lbf_variable` output.
Variant IDs are intersected in the first input's order. Defaults match
`coloc.susie`: p1 and p2 are 1e-4, p12 is 5e-6, and signal pairs with less
than 0.5 posterior overlap are omitted. `ColocOptions` can disable trimming,
change these priors, or provide positive per-variant prior weights.

The regression fixture is generated from coloc 5.2.3's bundled
`coloc_test_data` with susieR 0.14.2. Both the single-signal D1/D2 example and
the multi-signal D3/D4 example match R's H0-H4 and conditional-H4 variant
posteriors to numerical precision. Regenerate it with
`src/test/resources/r-reference/generate-coloc-susie-reference.R`.

## Observed-variable SEM from covariance

This path model estimates `x -> m -> y` and free residual variances:

```java
SemModel model = SemModel.builder("x", "m", "y")
    .regression("m", "x", 0.4)
    .regression("y", "m", 0.6)
    .variance("x", 1.0)
    .variance("m", 0.8)
    .variance("y", 0.5)
    .build();

double[] covariance = {
    1.0, 0.5, 0.35,
    0.5, 1.0, 0.70,
    0.35, 0.70, 1.0
};

SemFitResult sem = Sem.fitCovariance(
    covariance, 2_000, model,
    SemOptions.defaults(), BackendPolicy.PREFERRED);

SemParameterEstimate path = sem.parameter("y~m");
System.out.printf("beta=%g se=%g z=%g p=%g%n",
    path.estimate(), path.standardError(),
    path.zStatistic(), path.pValue());
```

Results include the likelihood chi-square, CFI, TLI, RMSEA, SRMR, AIC, BIC,
and implied covariance. A just-identified model can fit perfectly by
construction; fit indices are informative only with positive model degrees of
freedom.

## SEM from rows and equality constraints

```java
SemModel constrained = SemModel.builder("x", "m", "y")
    .regression("equalPath", "m", "x", 0.4)
    .regression("equalPath", "y", "m", 0.4)
    .build();

SemFitResult fromRows = Sem.fit(
    observedRows, constrained,
    SemOptions.defaults(), BackendPolicy.PREFERRED);
```

Repeated labels impose equality constraints. Use the builder's fixed
regression/variance/covariance methods for known parameters. The row-data fit
uses complete-case covariance ML according to `SemOptions`.

Current SEM scope is observed-variable covariance structure. Latent
measurement models, mean structures, ordinal thresholds, robust sandwich
corrections, modification indices, and FIML missingness are not silently
approximated. See the [scope document](../mr-timeseries-susie-sem.md).
