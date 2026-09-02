# Additive, additive mixed, and distributional models

JLinAlg compiles smooths into numeric basis and penalty blocks before fitting.
The iterative path therefore performs dense linear algebra through JDistlib
instead of interpreting formula terms. `BackendPolicy.PREFERRED` tries
auto-routed GPU execution first, then oneMKL, OpenBLAS, and portable Java CPU.

## Gaussian and generalized additive models

`Gam.fitGaussian` estimates ordinary P-spline smoothing parameters by exact
Gaussian REML. Select profile ML with
`RemlOptions.builder().varianceEstimation(VarianceEstimation.ML)`. The fast
default denominator DF is the documented residual approximation;
Satterthwaite and Kenward-Roger remain explicit options.

```java
PSplineTerm age = PSplineTerm.of("s(age)", ageYears, 12);
GamResult fit = Gam.fitGaussian(y, fixedDesign, List.of(age));

double beta = fit.parametricCoefficients()[1];
double se = fit.parametricStandardErrors()[1];
double t = fit.parametricTStatistics()[1];
double p = fit.parametricPValues()[1];
```

`GeneralizedGam.fit` uses penalized quasi-likelihood/IRLS with prior weights and
offsets for all `GlmFamilies`. `WeightedGam` supplies exact inverse-variance
weights and offsets for Gaussian REML.

For multiple penalties, `QuadraticSmoothTerm` is the common representation.
`GaussianSmoothSelector` supports fixed smoothing, GCV, UBRE, and AIC. Tensor
products retain one smoothing parameter per margin.

```java
QuadraticSmoothTerm surface = TensorProductPSplineTerm.of(
    "te(latitude,longitude)", latitude, longitude, 8, 8);
GaussianSmoothSelectionResult surfaceFit = GaussianSmoothSelector.fit(
    y, fixedDesign, List.of(surface));
double[] lambda = surfaceFit.smoothingParameters().get(0);
```

Available constructors include ordinary P-splines, cyclic Fourier/cubic
smooths, low-rank two-dimensional thin-plate smooths, tensor `te`, constrained
interaction `ti`, varying-coefficient terms, random-effect smooths, and
graph-Laplacian Markov random fields. `AdvancedGamFormula` recognizes `s`,
`te`, and `ti`:

```java
CompiledQuadraticGamFormula model = AdvancedGamFormula.compile(
    "expression ~ sex + te(age,bmi,kx=8,kz=7)", table);
GaussianSmoothSelectionResult result = model.fitGaussian(
    SmoothingSelectionOptions.gcv(), BackendPolicy.PREFERRED);
```

## High-throughput execution

`PreparedGam` retains univariate bases and warm-starts variance components for
repeated responses. `DiscretePSplineBasis` stores one row per unique covariate
value and forms `B'WB` and `B'Wy` without expanding repeated rows.

For GWAS/TWAS/EWAS/PWAS, `PreparedGamAssociation` retains the selected
penalized null-model factorization. A marker block requires one `X'G` product;
the remaining score/one-step effect calculations can use bounded CPU
parallelism. Results use the same effect-size, SE, statistic, regular p-value,
log10-p, and negative-log10-p contract as other association engines.

```java
try (PreparedGamAssociation scan = new PreparedGamAssociation(
        nullFit, AssociationEngineOptions.acceleratedSerial())) {
    AssociationBatchResult hits = scan.scan(genotypes, variantNames);
}
```

## GAMMs, pedigree/GRM covariance, and correlated errors

`Gamm.fitGaussian` combines smooth covariance blocks with any
`VarianceComponent`. `GammCovariances` directly aligns random intercepts,
random slopes, pedigrees, and genomic relationship matrices—including dense
cryptic-relatedness GRMs—to observation IDs. An optional residual correlation
matrix accepts AR/MA/ARMA structures from `ArmaCorrelation`.

`GeneralizedGamm` is the fast PQL compatibility route.
`GeneralizedGammLaplace` is the accuracy-oriented dense first-order Laplace
route. It eigendecomposes positive-semidefinite covariance bases once, profiles
their variances on the log scale, and returns conditional component modes plus
asymptotic beta/SE/z/p inference. Both accept pedigree and GRM components.

```java
VarianceComponent animal = GammCovariances.pedigree(
    "animal", observationIds, pedigree);
GeneralizedGammLaplaceResult fit = GeneralizedGammLaplace.fit(
    y, fixedDesign, List.of(ageSmooth), GlmFamilies.binomial(),
    List.of(animal), weights, offset,
    GlmmLaplaceOptions.defaults(), BackendPolicy.PREFERRED);
```

The Laplace implementation is a dense first-order approximation. PQL remains
useful for very large exploratory scans; exact Gaussian REML remains preferred
for Gaussian outcomes and finite-DF inference.

## GAMLSS- and VGAM-style models

`DistributionalFamily` declares any number of linked predictors. The block
Fisher solver estimates all parameters jointly, includes cross-parameter
information, step halving, quadratic smooth penalties, and JDistlib Cholesky
solves. Built-ins include Gaussian location/scale, Gamma mean/shape, Beta
mean/precision, NB2 mean/size, zero-inflated Poisson, hurdle Poisson,
multinomial logits, and adjacent-category ordinal logits.

```java
DistributionalResult fit = DistributionalFormula.fit(
    List.of("y ~ s(age,k=10)", "y ~ sex + s(age,k=7)"),
    List.of(List.of(new double[]{0.8}), List.of(new double[]{2.0})),
    table, DistributionalFamilies.gaussianLocationScale(),
    FormulaOptions.defaults(), DistributionalOptions.defaults(),
    BackendPolicy.PREFERRED);
```

For categorical vector responses, request normalized observation-by-category
probabilities from the same compiled parameter predictors:

```java
DistributionalFamily family = DistributionalFamilies.multinomial(3);
double[][] probability = DistributionalPrediction.categoryProbabilities(
    fit, parameterPredictors, family);
```

`PreparedDistributionalModel.refit` starts a nearby response at the preceding
coefficient blocks. `DistributionalPrediction` applies new compiled predictors;
`DistributionalDiagnostics` returns quantile residuals and Gaussian centiles;
`DistributionalModelComparison` supplies nested likelihood-ratio tests and AIC.

## Performance benchmark

Run `gradlew benchmarkGam` to time cold basis construction, compressed
cross-products, a fixed-penalty fit, a prepared marker scan, and a warm
distributional refit. The deterministic default workload uses 1,200
observations and 256 markers. A 2026-09-02 oneMKL acceptance run completed the
fixed fit in 2.214 seconds, scanned all markers in 0.773 seconds, and completed
the warm distributional refit in 0.730 seconds; the discrete basis retained
one-sixth as many rows and every scan result was valid. These are regression
observations, not cross-machine performance guarantees. Override workload size
with `-Djlinalg.benchmark.observations` and
`-Djlinalg.benchmark.markers`.

## R validation

Three deterministic, committed cross-language gates cover every stage without
requiring R during the Java build. `MgcvReferenceTest` compares Gaussian GAM
fitted values and EDF with `mgcv` 1.9-4. `Lme4LaplaceReferenceTest` compares a
binomial random-intercept Laplace fit with `lme4` 2.0-6.
`DistributionalRReferenceTest` compares Gaussian location/scale coefficients,
likelihood, and fitted parameters with `gamlss` 5.5-0, plus multinomial logits,
likelihood, and category probabilities with `VGAM` 1.1-14. The corresponding
generators and versioned properties are under
`src/test/resources/r-reference`. `generate-additive-optional-references.R`
also smoke-tests `gamm4` and `gam` when those packages are installed.
