# Repeated association, GWAS, and TWAS

## CLI-only association workflows

```bash
java -jar jlinalg-0.3.3.jar --pheno phenotype.csv --omics expression.csv --id SampleName --omics-type expression --formula "BMI ~ Sex + Age + <omics>" --transform "<omics> = winsor_mad(k=4) | zscore()" --annot genes.tsv --annot-id gene_id --annot-cols symbol,chromosome --out bmi-expression.csv

java -jar jlinalg-0.3.3.jar --pheno phenotype.csv --omics cohort.vcf.gz --id IID --formula "BMI ~ Sex + Age + PC1 + PC2 + <omics>" --min-maf 0.01 --min-mac 20 --max-marker-missing 0.02 --out bmi-gwas.csv

java -jar jlinalg-0.3.3.jar --pheno phenotype.csv --omics cohort.bgen --sample-file cohort.sample --id IID --formula "BMI ~ Sex + Age + PC1 + PC2 + <omics>" --grm cohort --variance-components null-model --out bmi-gwas-grm.tsv
```

BCF and biallelic layout-2 BGEN are supported. Numeric omics mixed scans refit
every feature; genotype LMM uses the explicitly logged null-model P3D/EMMAX
route. Add `--dry-run` to inspect alignment, model routing, block size, worker
capacity, and backend before scanning.

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

## Audited edge cases and output contracts

Predictor tests use scale-relative information checks, so changing a predictor's
units does not alone make it nonestimable. OLS computes alternative residual sums
of squares directly: subtracting two nearly equal sums can materially distort SEs
for strong signals. GLM/P3D estimates retain their documented score/null-model
estimands; they are not alternative-model refits. P3D rejects a nonconverged REML
null, and iterative built-in refit adapters propagate nonconvergence through the
selected failure policy instead of discarding it when extracting statistics.

VCF `DS` is usable even when no `GT` field is present. Alignment indices must be
nonnegative and unique, and a custom omics transform must preserve sample count.
Transforms receive a defensive copy. Cohort dosage variance uses centered updates
instead of cancellation-prone raw moments; NaN denotes missing genotype dosage,
whereas infinity is an error in variant QC.

Streaming sinks receive exactly one event per source item in source order.
`testedVariants`/`testedFeatures` count successful estimates, not failed attempts:
source = tested + excluded + failed (omics has no QC-excluded category). Collinear
features appear only as failures, not additional successful rows containing NaNs.

The CLI leaves unmatched annotation cells blank. Disk BH preserves quoted tabs,
quotes, and embedded newlines, uses at most 64 input runs in a merge pass, and
cleans scratch sort data on close. A `.partial` result remains recoverable after
an interrupted run; it is not a resumable checkpoint. BH adjusts over finite
p-values in [0,1]; failed/excluded rows do not enlarge that denominator. Define
the testing family before analysis rather than interpreting that convention as
a correction over every attempted/failed hypothesis.

Input/output/log aliases are rejected before the main CLI opens a log. Without
`--overwrite`, existing result/log/manifest/partial files are rejected. Existing
results cannot be authenticated by `--resume`, which now reports an error rather
than declaring any existing file complete. Use a new output path, or deliberately
rerun with `--overwrite` and without `--resume`. The standalone `mr-estimate`,
`beta-regression`, and `penalized-regression` commands require fresh output paths
(they do not expose an overwrite option).

For set tests, small positive weights no longer disappear behind an absolute
eigenvalue cutoff. Rank-two positive chi-square mixtures use an angular integral
over the Gaussian radius/angle decomposition. Factoring out the peak exponential
preserves relative accuracy for extreme tails and nearly rank-one kernels.
Simulation critical values invert the component tail rather than substituting
an unrelated moment-matched quantile. Zero padding does not change positive-rank
quantiles; near survival probability one, inversion evaluates the lower tail
directly. Higher-rank mixtures use a positive gamma expansion with a PGF-based
omitted-probability bound of 1e-12 relative to the result. An 8192-term resource
limit or unresolvable spectrum produces an explicit numerical error, not an
uncertified Imhof tail or silent moment fallback. This controls series truncation,
not all floating-point rounding. Analytic SKAT-O remains a moment-matched
calibration, not a claim of exact finite-sample calibration.
Underflowed mixture probabilities retain the existing `Double.MIN_VALUE` floor.

Fast OLS and P3D normalize predictors before imputation and cross-products and
rescale beta/SE afterward. This preserves inference at marker units as small as
1e-160 or as large as 1e160; it does not promise representable results for every
finite input. Direct OLS scans record unrepresentable exported estimates as
failures. P3D reports unavailable estimates as NaN.

The [audit evidence and runnable commands](../../src/benchmark/resources/pipeline-audit/AUDIT.md)
include R fixtures, edge-case tests, and raw warm timings. On the audited host,
400-sample/128-marker OLS took 0.751 ms versus R `lm.fit` at 5.20 ms, with maximum
beta/SE/p discrepancies of 2.5e-16/1.2e-16/1.7e-15. Disk BH took 78.11 ms for
20,000 output rows. R's in-memory BH took 1 ms,
but excludes Java's file I/O and is not an equivalent end-to-end speed comparison.
The rank-two integral took 0.0293 ms in Java and 0.0220 ms in R; accuracy, not a
speed win, motivates that repair. For rank three, lambda=(1,.5,.25), q=100,
positive gamma-series evaluation took 0.1147 ms in Java versus 3.20 ms in R;
both include coefficient generation and assert a relative remainder bound.
These are warm, workload-specific measurements. Reproduce with the registered
`benchmarkPipelineOlsAudit`, `benchmarkPipelineBhAudit`, and
`benchmarkPipelineKernelAudit` Gradle tasks and normal package tests, as linked.
