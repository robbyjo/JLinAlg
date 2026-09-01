# Cox proportional-hazards, frailty, and pedigree survival models

## Encode the survival response

Right-censored observations use one time and event flag per row:

```java
CoxSurvivalData survival = CoxSurvivalData.rightCensored(
    followUpTime, observedEvent);
```

For delayed entry or interval-form counting-process data, provide start, stop,
event, and optional nonnegative stratum IDs:

```java
CoxSurvivalData survival = new CoxSurvivalData(
    entryTime, exitTime, observedEvent, stratum);
```

Rows are at risk when `start < eventTime <= stop`. Every row must satisfy
`0 <= start < stop`; at least one event is required. Do not add an intercept
column to a Cox fixed-effect design because the baseline hazard absorbs it.

## Fit a fixed-effect Cox model

```java
CoxOptions options = CoxOptions.defaults()
    .withTies(CoxTies.EFRON);

CoxResult fit = CoxRegression.fit(
    survival, covariates, offset,
    options, BackendPolicy.PREFERRED);

double beta = fit.beta()[0];
double hazardRatio = fit.hazardRatios()[0];
double se = fit.standardErrors()[0];
double z = fit.zStatistics()[0];
double p = fit.pValues()[0];
double minusLogP = fit.negativeLog10PValues()[0];
```

The default Efron rule is generally preferable when tied event times are
present; Breslow is available explicitly. Results include ordinary and
`-log10` p-values, hazard-ratio confidence intervals, convergence diagnostics,
backend provenance, and stratum-specific baseline hazard/survival steps.

The right-censored fast path sorts once per stratum and accumulates risk-set
moments. Start-stop input uses a general reference path because membership can
both enter and leave as time changes.

## Add Gaussian shared frailty

Create grouped random intercepts using the existing sparse term builder, then
convert them to a Cox frailty term:

```java
CoxRandomEffectTerm center = CoxRandomEffectTerm.independent(
    RandomEffectTerm.randomIntercept("center", centerIds));

CoxMixedOptions mixedOptions = new CoxMixedOptions(
    CoxOptions.defaults(),
    new double[] {0.5}, // one initial variance per term
    30, 1e-4, 1e-8, 1e4);

CoxMixedResult mixed = CoxMixedModel.fit(
    survival, covariates, List.of(center), offset,
    mixedOptions, BackendPolicy.PREFERRED);

double centerVariance = mixed.randomEffects("center").variance();
double[] centerModes = mixed.randomEffects("center").modes();
```

Several independent or precision-defined Gaussian terms may be supplied. The
conditional modes use penalized partial likelihood. Variances are profiled by
a Laplace approximation that integrates the Gaussian frailties; this is not
REML. `fixef()` and `ranef()` aliases follow the rest of JLinAlg's mixed-model
surface.

## Adjust for cryptic relatedness

Use a caller-supplied or genotype-derived GRM as a named Gaussian frailty:

```java
GenomicRelationshipMatrix grm =
    GenomicRelationshipMatrix.fromVariants(
        relationshipVariants, sampleIds,
        GenomicRelationshipOptions.defaults(),
        BackendPolicy.PREFERRED);

CoxMixedResult relatedFit = CoxKinshipFrailty.fit(
    survival, covariates, observationSampleIds, grm,
    offset, mixedOptions, 1e-8, BackendPolicy.PREFERRED);

double variance = relatedFit.randomEffects("kinship").variance();
double sampleMode = relatedFit.randomEffects("kinship").mode("sample-1");
```

Observation IDs may repeat for start-stop rows. The explicit diagonal
regularization supports empirical GRMs that are singular because of duplicate
samples or finite marker rank. The current covariance inversion and Laplace
random-information block are dense.

## Use pedigree-correlated frailty

Observation IDs may repeat, and pedigree individuals without survival rows
remain in the Gaussian frailty system:

```java
CoxPedigreeResult pedigreeFit = CoxPedigreeFrailty.fit(
    survival, covariates, observationAnimalIds, pedigree,
    offset, mixedOptions, BackendPolicy.PREFERRED);

double additiveFrailtyVariance = pedigreeFit.frailtyVariance();
double offspringMode = pedigreeFit.frailty("offspring");
Map<String, Double> allModes = pedigreeFit.ranef();
```

The model uses the pedigree's directly constructed sparse `A^-1` as its unit
Gaussian precision, avoiding a numerical inverse of the relationship matrix.
The current Laplace solve materializes the random-effect information block, so
very large pedigrees still need the planned sparse Cox precision engine.

## Interpretation and current boundaries

- `exp(beta)` is a covariate hazard ratio under proportional hazards.
- Fixed-effect p-values are asymptotic Wald z tests.
- Baseline survival is `exp(-cumulativeHazard)` within each stratum.
- Mixed, GRM, and pedigree estimates use Gaussian log frailty and a Laplace
  approximation. Gamma frailty and adaptive quadrature are not implemented.
- Schoenfeld residual proportional-hazards tests, martingale/deviance
  residuals, robust cluster sandwich covariance, recurrent-event robust
  inference, time-varying coefficient builders, and high-throughput prepared
  Cox score scans remain future work.

Numerical regression tests lock fixed Efron, Breslow, and delayed-entry results
to an independent `statsmodels` Cox implementation. Mixed and pedigree tests
exercise positive variance profiles, hazard-ratio inference, named modes, and
a non-identity parent-offspring relationship.
