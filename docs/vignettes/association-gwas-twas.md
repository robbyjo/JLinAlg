# Repeated association, GWAS, and TWAS

JLinAlg distinguishes fast prepared scans from exact repeated fits. Both
return ordered effect size, SE, statistic, p-value, log10 p-value, and
`-log10(p)`.

For production file-backed work, use `VariantSources` with CSV, TSV, VCF,
VCF.gz, BCF, or BGEN. PLINK is intentionally outside this API. Sample alignment
precedes cohort-specific MAC/MAF/missingness/quality filtering.

```java
VariantSource source = VariantSources.open(Path.of("cohort.bgen"));
VariantFilterOptions qc = VariantFilterOptions.builder()
    .minimumMac(20).minimumMaf(0.005)
    .maximumMissingRate(0.02)
    .minimumImputationQuality(0.8).build();

try (DelimitedAssociationWriter writer =
        new DelimitedAssociationWriter(Path.of("gwas.tsv"), '\t')) {
    StreamingAssociationPipeline.fastOlsTo(
        source, analysisSampleIds, phenotype, covariates,
        null, null, OlsOptions.defaults(), execution,
        new AssociationPipelineOptions(4096, qc), writer);
}
```

The sink overload keeps input and output memory bounded and writes results,
QC exclusions, and failures in source order.

## Fast OLS predictor scan

Rows are samples, covariate columns are shared by every test, and candidate
columns are markers or molecular features:

```java
AssociationEngineOptions execution = AssociationEngineOptions.cpuParallel()
    .withParallelism(8)
    .withChunkSize(256);

AssociationBatchResult scan = FastOlsAssociation.scanPredictors(
    phenotype,
    covariates,
    dosages,
    markerNames,
    null,                 // weights
    null,                 // offset
    OlsOptions.defaults(),
    execution);

AssociationEstimate hit = scan.estimate(0);
System.out.printf("%s beta=%g se=%g p=%g -log10p=%g%n",
    hit.name(), hit.beta(), hit.standardError(), hit.pValue(),
    hit.negativeLog10PValue());
```

`FastOlsAssociation` factors the covariates once and uses
Frisch-Waugh-Lovell residualization for marker blocks. The default changing-
predictor missing policy mean-imputes each marker. Select
`VariableMissingPolicy.ERROR` when imputation is not scientifically intended.

To scan many phenotypes against one design, call `scanResponses` and specify
the fixed-effect coefficient index to report.

## Fast GLM score scan

Prepare one null IRLS model, then score many predictors:

```java
FastGlmAssociation prepared = FastGlmAssociation.prepare(
    binaryPhenotype, covariates, GlmFamilies.binomial(),
    null, null, GlmOptions.defaults(), execution);

AssociationBatchResult scoreScan = prepared.scan(
    dosages, markerNames, execution);
```

The returned beta is a one-step null-model score estimate. Use it for screening
and exactly refit selected variants if the final reported beta must be the
full alternative-model maximum-likelihood estimate.

## Exact parallel refits

The generic engine appends each candidate as the last fixed-effect column and
fits the selected adapter independently:

```java
AssociationFitter fitter = AssociationModels.reml(
    varianceComponents, RemlOptions.defaults());

AssociationBatchResult exact = ParallelAssociationEngine.scanPredictors(
    phenotype, covariates, dosages, markerNames, fitter, execution);
```

Built-in adapters cover OLS/weighted OLS, GLM, dense or sparse LMM, correlated
LMM, REML, pedigree REML, GLMM PQL, pedigree GLMM PQL, ridge, and ARIMA-error
LMM. Implement `AssociationFitter` for another thread-safe model returning
`AssociationStatistics`.

For a long production scan, consider:

```java
AssociationEngineOptions tolerant = execution.withFailurePolicy(
    AssociationFailurePolicy.RECORD_NAN);
```

This preserves structured failures and input order rather than terminating the
whole job. Inspect `failures()` before writing final results.

## P3D/EMMAX-style mixed-model scan

When frozen null-model variance components are acceptable, prepare the REML
projection once:

```java
GenomicRelationshipMatrix grm =
    GenomicRelationshipMatrix.fromVariants(
        relationshipVariants, sampleIds,
        GenomicRelationshipOptions.defaults(),
        BackendPolicy.PREFERRED);

RemlAssociationScanner prepared = RemlAssociationScanner.prepare(
    phenotype, covariates, List.of(
        grm.varianceComponent("cryptic"),
        VarianceComponent.identity("residual", sampleIds.size())),
    RemlOptions.defaults(), BackendPolicy.PREFERRED);

AssociationScanResult p3d = prepared.scan(
    rowMajorDosages, markerCount, markerNames,
    new AssociationScanOptions(
        2048, GenotypeMissingPolicy.MEAN_IMPUTE, 1));
```

This does not re-estimate variance components for every marker. It is therefore
much faster than exact REML refits but represents the P3D/EMMAX approximation.
Increase scan parallelism primarily for a CPU backend; with GPU or native BLAS,
start with one submitting thread.

## Output hygiene

- Preserve allele/effect-direction metadata outside the numerical matrix.
- Confirm sample order is identical across phenotype, covariates, kinship, and
  candidate columns.
- Report the tested coefficient, model/approximation, missing policy, and
  backend provenance.
- Use `negativeLog10PValue()` for plotting but retain ordinary p-values and
  effect/SE for downstream meta-analysis.

See the [association engine reference](../association-engine.md) for exact
adapter signatures and execution-policy details.

## Omics transforms and rare-variant sets

`DelimitedMatrixSource` and `StreamingOmicsAssociationPipeline` provide the
same sample alignment for TWAS, EWAS, and PWAS feature matrices. Compose
Winsorization, log, z-score, or tie-aware rank inverse-normal transforms before
the selected missing-data policy.

Burden, SKAT, and SKAT-O accept explicit weighted `VariantSet` membership.
`LinearSetTestNullModel` handles unrelated samples; `RemlSetTestNullModel`
reuses the fitted GRM-adjusted mixed projection for related samples, including
SKAT-O calibration. The complete file and set-test contract is in the
[pipeline guide](../gwas-twas-pipeline.md).
