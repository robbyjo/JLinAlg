# GWAS, TWAS, EWAS, and PWAS pipeline

This layer joins file-backed inputs, aligned-cohort QC, transformations, and
prepared association engines. It deliberately does not accept PLINK files.
Genotypes remain block-streamed and `fastOlsTo`/`remlP3dTo` can stream results
directly to a sink, so neither variants nor output rows need to accumulate in
memory.

## Inputs

`VariantSources.open(path)` infers these suffixes:

| Input | Current contract |
| --- | --- |
| CSV / TSV | one variant per row; `id` followed by samples, or `id,chromosome,position,ref,alt` followed by samples |
| VCF / VCF.gz | `DS` dosage preferred, `GT` fallback; one row per alternate allele; INFO/R2/DR2 quality recognized |
| BCF | BCF 2.1/2.2 through HTSJDK-compatible codecs |
| BGEN | layout 2, biallelic phased or unphased, variable ploidy, 1--32 bit probabilities, uncompressed/zlib/zstd |

BGEN files without an embedded sample block require IDs explicitly:

```java
VariantSource variants = VariantSources.open(
    Path.of("cohort.bgen"), DataFormat.BGEN, analysisSampleIds);
```

Multiallelic VCF/BCF records expand by alternate allele. Multiallelic BGEN and
legacy BGEN layout 1 currently fail explicitly rather than being decoded
incorrectly. Random-access indexes and region queries are not implemented yet.

## Aligned-cohort QC

Sample IDs are matched exactly and reordered before statistics are calculated.
MAF, MAC, missingness, and imputation-quality filters therefore describe the
actual analysis cohort, not an unaligned source file. Missing values are
counted before any scan-time mean imputation.

```java
VariantFilterOptions qc = VariantFilterOptions.builder()
    .minimumMac(20)
    .minimumMaf(0.005)
    .maximumMissingRate(0.02)
    .minimumImputationQuality(0.8)
    .build();

AssociationPipelineOptions pipeline =
    new AssociationPipelineOptions(4096, qc);
```

There is no implicit rare-variant cutoff. Defaults exclude only monomorphic
variants; every other threshold must be selected by the analysis.

## Bounded-memory GWAS

For an unrelated quantitative trait, the fastest implemented path factors the
shared covariates once and residualizes genotype blocks:

```java
VariantSource source = VariantSources.open(Path.of("cohort.vcf.gz"));
AssociationEngineOptions execution =
    AssociationEngineOptions.acceleratedSerial()
        .withBackendPolicy(BackendPolicy.PREFERRED);

try (DelimitedAssociationWriter output =
        new DelimitedAssociationWriter(Path.of("gwas.tsv"), '\t')) {
    AssociationPipelineSummary summary =
        StreamingAssociationPipeline.fastOlsTo(
            source, analysisSampleIds, phenotype, covariates,
            weights, offset, OlsOptions.defaults(), execution,
            pipeline, output);
}
```

The output includes variant metadata, called/missing counts, MAF/MAC, beta,
SE, t statistic, denominator DF, p-value, `log10(p)`, `-log10(p)`, exclusions,
and structured failures. Use `fastOls` instead when an in-memory result is
convenient for a smaller analysis.

For cryptic relatedness, build a GRM from aligned, preferably LD-pruned
variants, inspect related pairs, and reuse it as a variance component:

```java
GenomicRelationshipMatrix grm =
    GenomicRelationshipMatrix.fromVariants(
        relationshipVariants, analysisSampleIds,
        GenomicRelationshipOptions.defaults(),
        BackendPolicy.PREFERRED);

List<RelatednessPair> relatives = grm.relatedPairsByKinship(0.0884);
List<VarianceComponent> kinshipComponents = List.of(
    grm.varianceComponent("grm"),
    VarianceComponent.identity("residual", grm.samples()));
```

GRM construction mean-imputes missing dosages after call-rate/MAF filtering
and routes the standardized cross-product through JDistlib. For repeated
phenotype rows, use `grm.varianceComponent("grm", observationSampleIds)` to
expand the sample covariance to observation scale.

Then prepare variance components once and stream a
P3D/EMMAX-style scan:

```java
RemlAssociationScanner nullModel = RemlAssociationScanner.prepare(
    phenotype, covariates, kinshipComponents,
    RemlOptions.defaults(), BackendPolicy.PREFERRED);

try (DelimitedAssociationWriter output =
        new DelimitedAssociationWriter(Path.of("mixed-gwas.tsv"), '\t')) {
    StreamingAssociationPipeline.remlP3dTo(
        source, analysisSampleIds, nullModel,
        new AssociationScanOptions(
            4096, GenotypeMissingPolicy.MEAN_IMPUTE, 1),
        pipeline, output);
}
```

This reuses frozen null-model variance estimates. Exact per-variant REML fits
remain available through `ParallelAssociationEngine` but cost much more.

## TWAS, EWAS, and PWAS matrices

`DelimitedMatrixSource` reads a CSV/TSV feature-by-sample matrix in blocks.
`StreamingOmicsAssociationPipeline.scanResponses` scans changing molecular
traits, while `scanPredictors` scans changing molecular predictors. Available
transforms preserve missing values until the selected missing-data policy:

```java
OmicsTransform transform = OmicsTransforms.compose(
    OmicsTransforms.winsorize(0.01, 0.99),
    OmicsTransforms.rankInverseNormal());

OmicsAssociationResult results =
    StreamingOmicsAssociationPipeline.scanResponses(
        expression, analysisSampleIds, design, testedCoefficient,
        transform, OmicsMissingPolicy.MEAN_IMPUTE, 512,
        weights, offset, OlsOptions.defaults(), execution);
```

Identity, empirical-quantile Winsorization, `log1p`, shifted log, sample
z-score, tie-aware Blom rank inverse-normal, and transform composition are
implemented. These matrix scans currently return in-memory results; an
incremental omics sink is still planned.

## Burden, SKAT, and SKAT-O

Build a `VariantSet` from aligned `WeightedVariant` members. Effect-allele
orientation, missingness policy, and the same MAF/MAC filters are applied
before testing. Conventional Beta(1,25) weights are available through
`VariantWeights.betaBurden` and `betaKernel`.

```java
LinearSetTestNullModel nullModel = LinearSetTestNullModel.prepare(
    phenotype, covariates, OlsOptions.defaults(), BackendPolicy.PREFERRED);

SetTestResult burden = SetTests.burden(
    gene, nullModel, SetTestOptions.defaults());
SetTestResult skat = SetTests.skat(
    gene, nullModel, SetTestOptions.defaults());
SkatOResult skatO = SetTests.skatO(
    gene, nullModel, SetTestOptions.defaults());
```

SKAT uses a positive chi-square mixture with deterministic Imhof integration
and guarded moment matching. SKAT-O reports every rho-grid component and a
seeded correlated-null adjusted p-value. For related samples, wrap the retained
mixed projection once:

```java
RemlSetTestNullModel related = new RemlSetTestNullModel(nullModelScanner);
SetTestResult mixedBurden = SetTests.burden(gene, related, options);
SetTestResult mixedSkat = SetTests.skat(gene, related, options);
SkatOResult mixedSkatO = SetTests.skatO(gene, related, options);
```

Gene/region annotation, window construction, conditional analysis, and
external reference-panel LD are the next orchestration layer; callers currently
provide set membership explicitly.

## Current boundaries

- Quantitative-trait single-variant streaming is implemented for fast OLS and
  frozen-null REML. File-backed binary/count GLM streaming still needs wiring
  to the prepared GLM score engine.
- BGEN multiallelic decoding, `.bgi`/tabix/CSI region queries, and resume
  checkpoints remain open.
- PQL is the current fast GLMM approximation. Exact `glmer`-class Laplace and
  adaptive Gauss-Hermite likelihoods are not implemented.
- Full `lme4`/`pedigreemm` API and likelihood parity is not yet claimed; see
  [the compatibility roadmap](lme4-pedigreemm-roadmap.md).

The source design was informed by the user-authorized GPU_eQTL repository, but
this implementation is clean-room rather than copied source so JLinAlg can
retain its GPL-2.0-or-later licensing.
