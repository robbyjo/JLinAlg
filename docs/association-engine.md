# Repeated association engine

JLinAlg has ordered, bounded-parallel scans for the two common repeated-fit
shapes:

- one response and many changing predictor columns;
- many response columns and one fixed design.

Every result contains the input name, beta, standard error, t or z statistic,
p-value, and denominator degrees of freedom. `beta()` and `effectSizes()` are
equivalent. Results remain in input order even when fits complete out of order.

Regular, literal `log10(p)`, and the GWAS-standard `-log10(p)` are all retained:

```java
double[] p = scan.pValues();
double[] logP = scan.pValues(PValueScale.LOG10);
double[] minusLogP = scan.pValues(PValueScale.NEGATIVE_LOG10);
```

The logarithmic forms are evaluated directly in the distribution's log-tail,
so they remain useful when a regular floating-point p-value underflows to zero.

## Fast OLS

`FastOlsAssociation` factors the shared covariate design once. Predictor scans
then residualize blocks of predictors with Frisch-Waugh-Lovell algebra, while
response scans reuse the fixed `X'X` inverse. Both paths support positive
weights and an additive offset.

```java
AssociationEngineOptions execution = AssociationEngineOptions.cpuParallel()
    .withParallelism(8)
    .withChunkSize(256);

AssociationBatchResult scan = FastOlsAssociation.scanPredictors(
    phenotype,
    covariates,          // observation by fixed-covariate matrix
    dosages,             // observation by marker matrix
    markerNames,
    weights,             // null means all one
    offset,              // null means all zero
    OlsOptions.defaults(),
    execution);

AssociationEstimate first = scan.estimate(0);
double beta = first.beta();
double se = first.standardError();
double t = first.statistic();
double p = first.pValue();
double minusLog10P = first.negativeLog10PValue();
```

Missing changing predictors can be mean-imputed or rejected. The fixed design,
response, weights, and offset must be complete on the fast path. Use the generic
engine with the model's missing-data policy when trait-specific complete-case
fitting is required.

## Exact repeated fits with any model

`ParallelAssociationEngine` performs an independent model fit for every
changing predictor or response. A predictor is appended as the final design
column, so its coefficient index is the number of shared fixed columns.

```java
AssociationFitter model = AssociationModels.reml(
    varianceComponents, remlOptions);

AssociationBatchResult scan = ParallelAssociationEngine.scanPredictors(
    phenotype, covariates, dosages, markerNames, model, execution);
```

Changing-response scans select the fixed-effect coefficient to return:

```java
AssociationBatchResult traits = ParallelAssociationEngine.scanResponses(
    expressionTraits, design, testedCoefficient, traitNames,
    AssociationModels.ols(OlsOptions.defaults()), execution);
```

Built-in thread-safe adapters include:

- OLS and weighted OLS;
- GLMs;
- exact REML and Gaussian LMMs, including a known residual correlation;
- pedigree REML;
- PQL GLMMs, with or without pedigree covariance;
- ridge regression with effective-DF model-based inference;
- profile-REML LMMs with ARIMA errors.

Any other model can participate by implementing the small `AssociationFitter`
contract and returning `AssociationStatistics`. A fitter must be thread-safe
and must not mutate or retain the supplied response and design work arrays.

## Prepared-null GLM score scan

`FastGlmAssociation` fits the shared-covariate null GLM once and caches the
final IRLS working weights, fixed-effect information inverse, and score
residual. Predictor blocks are residualized in that metric and return an
efficient-score z test, one-step effect estimate, SE, and regular/log p-values.

```java
FastGlmAssociation prepared = FastGlmAssociation.prepare(
    phenotype, covariates, GlmFamilies.binomial(),
    priorWeights, offset, GlmOptions.defaults(), execution);

AssociationBatchResult scan = prepared.scan(
    dosages, markerNames, execution);
```

The one-step beta is a null-model score estimate, not an exact alternative-model
maximum-likelihood coefficient. Use it for high-throughput screening and run
exact GLM refits when the final reported effect must be the full alternative
model estimate.

## Fast REML/P3D scan

For GWAS/TWAS where frozen null-model variance components are acceptable,
`RemlAssociationScanner` fits REML once, caches the GLS projection, and scans
marker blocks. This is the P3D/EMMAX-style path and is substantially cheaper
than exact per-marker REML refits.

```java
RemlAssociationScanner prepared = RemlAssociationScanner.prepare(
    phenotype, covariates, varianceComponents, remlOptions,
    BackendPolicy.PREFERRED);

AssociationScanResult scan = prepared.scan(
    rowMajorDosages, markerCount, markerNames,
    new AssociationScanOptions(
        2048, GenotypeMissingPolicy.MEAN_IMPUTE, 1));
```

Set the last option above to the desired batch parallelism for a CPU backend.

## Parallel and accelerator policy

`AssociationEngineOptions.defaults()` (also named
`acceleratedSerial()`) uses `BackendPolicy.PREFERRED` and one submitting thread,
allowing JDistlib to route large matrix operations to GPU, then oneMKL,
OpenBLAS, and finally portable CPU. This is the default because it follows the
project-wide accelerator preference without multiplying Java workers by
BLAS/GPU workers.

`AssociationEngineOptions.cpuParallel()` instead selects portable CPU and one
outer Java worker per available processor. Use `withParallelism()` to cap it.

`chunkSize` bounds per-worker setup and scheduling overhead. Larger chunks are
usually better for inexpensive fits; smaller chunks improve load balancing for
mixed models whose convergence times differ.

`FAIL_FAST` aborts on the first failed fit. `RECORD_NAN` retains a structured
failure and leaves that result's numeric fields as `NaN`, allowing long scans
to finish. The engine never silently reorders or drops inputs.

LASSO and elastic-net coefficient p-values are intentionally not exposed as a
generic association adapter: ordinary post-selection p-values are not valid
without an explicitly chosen inferential procedure. Ridge has a built-in
model-based adapter; callers can supply a custom selection-adjusted fitter for
other penalized analyses.
