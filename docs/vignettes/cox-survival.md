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
moments. Start-stop input builds one counting-process plan and sweeps event
times in ascending order, adding rows with `start < time` and removing rows
with `stop < time`. Both paths are reused by `CoxRegression.Prepared`.

## Robust inference and diagnostics

```java
CoxDiagnosticsResult diagnostics = CoxDiagnostics.analyze(
    survival, covariates, fit, offset, recurrentSubjectIds);

double[] robustSe = diagnostics.robustStandardErrors();
CoxProportionalHazardsTest ph = diagnostics.proportionalHazardsTest();
CoxResiduals residuals = diagnostics.residuals();
```

The sandwich meat is aggregated by cluster, supporting repeated start-stop
rows. Residual exports include martingale, deviance, per-row score, dfbeta,
and event-row Schoenfeld residuals. The proportional-hazards interface returns
per-term and global score tests using log event time.

## Prepared score scans

`FastCoxAssociation.prepare` fits the null model once and reuses its risk-set
plan for blockwise model-based, cluster-robust, or caller-correlation-aware
score tests. `StreamingOmicsAssociationPipeline.scanPredictorsCoxTo` connects
the same engine to file-backed blocks and an incremental result sink. Run
`benchmarkCoxPipeline` for prepared-score/full-refit and start-stop timings.
On the 2026-09-05 development-host run (2,000 rows, 512 predictors), the
prepared engine processed 28,748 predictors/second versus 1,829 for
single-threaded R `survival::coxph.fit` refits, a 15.72x throughput advantage.
The R comparator is `src/benchmark/r/cox_pipeline_benchmark.R`.

## Shared gamma frailty

```java
CoxGammaFrailtyResult gamma = CoxGammaFrailty.fitAtVariance(
    survival, covariates, centerIds, offset, 0.2,
    CoxGammaFrailtyOptions.defaults(), BackendPolicy.PREFERRED);
```

Frailties are multiplicative gamma variables with mean one and variance
`theta`; logged conditional modes enter the linear predictor. Fixed-`theta`
coefficients, SEs, and modes are regression-tested against
`survival::coxph(... + frailty(..., distribution="gamma", method="fixed"))`.
`fit` profiles `theta` with a Laplace approximation on its log scale; that
explicit approximation need not equal R's gamma-specific EM/profile criterion.

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
The default `fit` path materializes the random-effect information block.

For one-stratum right-censored data with distinct event times, opt into the
sparse precision solve:

```java
PedigreeRandomEffectTerm sparseTerm =
    PedigreeRandomEffectTerm.ofSparse(
        "pedigree", observationAnimalIds, pedigreeEntries,
        inbreedingCoefficients);

CoxPedigreeResult sparsePedigree = CoxPedigreeFrailty.fitSparse(
    survival, covariates, sparseTerm,
    offset, mixedOptions, BackendPolicy.PREFERRED);

CoxMixedResult diagnostics = sparsePedigree.mixedModel();
int equationNnz = diagnostics.sparseEquationNonzeroCount();
int factorNnz = diagnostics.sparseFactorNonzeroCount();
boolean boundaryFit = diagnostics.isSingular(1e-6);
```

Building the term with `ofSparse` (or `ofUninbred` when justified) avoids
materializing dense `A`; the convenience `fitSparse` overload accepting an
already-built `Pedigree` only makes the Cox equation solve sparse. The sparse
path supports repeated rows per pedigree individual and retains unobserved
ancestors. `SparseCoxMixedModel` can also consume a generic
unit-incidence `RandomEffectTerm` and caller-supplied `SparsePrecisionMatrix`,
which covers independent shared frailty and externally prepared sparse
precision models. It keeps the exact penalized score, but approximates the
profiled random-effect information by its diagonal plus the full sparse
precision. Consequently conditional modes agree with the dense score solution
when both converge, while Laplace-profiled variances and fixed-effect covariance
need not be identical. Use the dense path when the approximation or sparse
kernel restrictions are unsuitable.

## Interpretation and current boundaries

- `exp(beta)` is a covariate hazard ratio under proportional hazards.
- Fixed-effect p-values are asymptotic Wald z tests.
- Baseline survival is `exp(-cumulativeHazard)` within each stratum.
- Mixed, GRM, and pedigree estimates use Gaussian log frailty and a Laplace
  approximation. Gamma frailty uses its separately documented multiplicative
  gamma/Laplace path. Adaptive quadrature is not implemented.
- `solver()` distinguishes dense and sparse-precision results. Sparse results
  additionally report coefficient, equation-nonzero, and factor-nonzero counts;
  `isSingular(tolerance)`, `converged()`, `convergenceMessage()`, and `backend()`
  expose the boundary, optimization, and compute-backend state.
- Time-varying coefficient builders remain caller-defined; diagnostics expose
  the evidence needed to motivate them.

Numerical regression tests lock fixed Efron, Breslow, and delayed-entry results
to an independent `statsmodels` Cox implementation. Mixed and pedigree tests
exercise positive variance profiles, hazard-ratio inference, named modes, a
non-identity parent-offspring relationship, repeated sparse incidence, and
dense/sparse conditional-mode parity.
