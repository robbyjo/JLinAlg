# Parallel xWAS Mendelian-randomization pipeline

High-throughput molecular MR often has two independent dimensions: many
exposures such as genes, CpGs, proteins, or metabolites, and many outcome
phenotypes organized into disease families. `XwasMrPipeline` evaluates that
complete exposure-by-outcome grid while reusing each exposure's already-
clumped instruments.

## Why the pipeline is staged

A direct nested loop calling `MendelianRandomization.analyze` in every cell is
wasteful. The combined analysis runs a weighted-median bootstrap and
leave-one-out fits even for clearly null pairs. It can also create millions of
queued futures and make output order depend on thread scheduling.

The prepared pipeline instead:

1. accepts immutable, already-clumped instruments for every exposure;
2. indexes every outcome once by variant ID;
3. looks up only the exposure instrument SNPs for each phenotype;
4. harmonizes each exposure-outcome pair;
5. runs fast fixed or multiplicative-random IVW screening;
6. runs the complete MR diagnostics only for pairs passing the threshold;
7. processes pairs in bounded parallel blocks; and
8. merges hits in exposure-major, outcome-major order regardless of thread
   count.

Clumping is intentionally outside the nested loop. If an exposure's candidate
set, LD database, ancestry panel, and clumping settings have not changed, its
clumped instrument set is reused for every outcome.

## Prepare exposures once

Use the instrument search/format and clump commands described in the
[end-to-end MR vignette](mr-end-to-end.md). Convert each clumped exposure table
to `SummaryAssociation` objects and retain its molecular-feature identifier:

```java
List<XwasMrExposure> exposures = List.of(
    new XwasMrExposure(
        "ENSG00000123456",
        "GENE1 expression, whole blood",
        gene1ClumpedAssociations),
    new XwasMrExposure(
        "cg01234567",
        "cg01234567 methylation",
        cpgClumpedAssociations),
    new XwasMrExposure(
        "P12345",
        "Protein P12345",
        proteinClumpedAssociations));
```

`clumpedInstruments` contains exposure beta, SE, effect allele, other allele,
and preferably effect-allele frequency. Do not clump again for each outcome.
Create a new `XwasMrExposure` only when the exposure instrument definition
changes.

## Organize many phenotypes

An outcome carries both a unique ID and an optional category. Categories are
metadata for reporting and stratification; all outcomes participate in the
same efficient scan:

```java
List<XwasMrOutcome> outcomes = List.of(
    new XwasMrOutcome(
        "CAD", "Coronary artery disease", "cardiovascular",
        cadAssociations),
    new XwasMrOutcome(
        "LDL", "LDL cholesterol", "cardiovascular-risk-factor",
        ldlAssociations),
    new XwasMrOutcome(
        "CKD", "Chronic kidney disease", "kidney",
        ckdAssociations),
    new XwasMrOutcome(
        "eGFR", "Estimated glomerular filtration rate",
        "kidney-risk-factor", egfrAssociations),
    new XwasMrOutcome(
        "COPD", "Chronic obstructive pulmonary disease", "lung",
        copdAssociations),
    new XwasMrOutcome(
        "FEV1", "Forced expiratory volume in one second",
        "lung-risk-factor", fev1Associations));
```

Each outcome association list may contain a genome-wide table. Preparation
indexes it once, after which every pair performs lookups proportional to the
number of clumped instruments rather than rescanning the whole outcome.
Outcome IDs must be unique even when labels repeat across cohorts.

Phenotype selection remains a scientific design decision. Define disease and
risk-factor families before looking at MR results, document exclusions and
sample overlap, and avoid choosing phenotypes only because preliminary results
are significant.

## Run the complete scan from the CLI

The `mr-xwas` command exposes the same bounded, parallel pipeline without
requiring Java code. Supply one long-format table of already-clumped exposure
instruments and one long-format outcome table containing all phenotypes:

```console
java -jar jlinalg-<version>.jar mr-xwas --exposure molecular-traits.clumped.tsv --outcome phenotypes.tsv.gz --output xwas-mr-results.tsv --fdr-output xwas-mr-all-pairs.tsv --p-threshold 5e-8 --threads 8 --pair-block-size 256
```

Each exposure row needs a SNP, beta, SE, effect allele, other allele, and an
exposure identifier. Each outcome row needs the corresponding association
columns and an outcome identifier. The CLI recognizes common TwoSampleMR-style
names such as `SNP`, `beta.exposure`, `se.exposure`, `effect_allele.exposure`,
`id.exposure`, and their `.outcome` counterparts. Use
`--exposure-id-column` or `--outcome-id-column` when identifiers use custom
names. An optional `category` or `phenotype_category` outcome column preserves
disease families such as cardiovascular, kidney, and lung in the result.

Exactly one screening threshold is required:

- `--p-threshold 5e-8` retains `p <= 5e-8`;
- `--log10-p-threshold -7.30103` retains `log10(p) <= -7.30103`;
- `--negative-log10-p-threshold 7.30103` retains `-log10(p) >= 7.30103`.

Inputs and outputs may be CSV, TSV, or gzip-compressed. Multiplicative-random
IVW is the default screen; add `--screen-method fixed` for fixed-effect IVW.
The output contains the screening estimate and heterogeneity statistics plus
hit-only MR-Egger, weighted-median, instrument-strength, harmonization, and
warning fields. A sibling failures table is written automatically; use
`--failures FILE` to choose its location. Run `mr-xwas --help` for all options,
including bootstrap replicates, confidence level, seed, and overwrite control.

When `--fdr-output FILE` is supplied, the CLI also writes an uncompressed TSV
containing every successfully screened pair in deterministic grid order. Its
`threshold_passed` field identifies pairs sent to full diagnostics, and its
`fdr_bh` field is adjusted over all finite screening p-values in that file.
Pairs with insufficient instruments cannot be screened and are not members of
the BH family; screening or diagnostic exceptions remain explicit in the
failure table.

See [xWAS MR CLI and performance](../xwas-mr-cli-performance.md) for the full
input/output contract, reproducible benchmark commands, numerical comparison
with R, and measured throughput.

## Make the threshold scale unambiguous

Three mathematically equivalent threshold forms have different comparison
directions:

```java
XwasMrSignificanceFilter byP =
    XwasMrSignificanceFilter.pValueAtMost(5e-8);          // p <= 5e-8

XwasMrSignificanceFilter byLogP =
    XwasMrSignificanceFilter.log10PAtMost(Math.log10(5e-8));
                                                        // log10(p) <= -7.301...

XwasMrSignificanceFilter byMinusLogP =
    XwasMrSignificanceFilter.negativeLog10PAtLeast(
        -Math.log10(5e-8));                             // -log10(p) >= 7.301...
```

Thus “log10P less than a threshold” means `log10(p) <= threshold`. For the
more common positive Manhattan-plot scale, significance means `-log10(p) >=
threshold`, not less than it. The filter constructors reject thresholds with
an impossible sign.

The threshold should reflect the number and dependency of the planned tests.
For a simple Bonferroni family-wise threshold over the complete grid:

```java
long tests = Math.multiplyExact((long) exposures.size(), outcomes.size());
double alpha = 0.05 / tests;
XwasMrSignificanceFilter filter =
    XwasMrSignificanceFilter.pValueAtMost(alpha);
```

For FDR, pass `--fdr-output xwas-mr-all-pairs.tsv`. The disk-backed BH pass uses
every successfully screened pair, including pairs that do not pass the primary
diagnostic threshold, and appends `fdr_bh` without changing bounded scan memory.
Filtering first at nominal 0.05 and then applying BH only to retained rows would
be invalid; the all-pairs file avoids that error.

## Run the two-dimensional scan

```java
XwasMrPipeline pipeline = XwasMrPipeline.prepare(exposures, outcomes);

int threads = Math.min(16,
    Math.max(1, Runtime.getRuntime().availableProcessors()));

XwasMrOptions options = new XwasMrOptions(
    threads,
    256, // maximum exposure-outcome pairs resident in one work block
    XwasMrScreeningMethod.IVW_MULTIPLICATIVE_RANDOM,
    filter,
    new MrOptions(0.95, 1_000, 20260904L));

XwasMrBatchResult result = pipeline.scan(options);
```

`pairBlockSize` controls temporary result memory, not statistical batching.
A value between roughly `4 * threads` and `32 * threads` is a reasonable
starting point. Increase it only after profiling. The number of worker threads
should normally be bounded by available processors and memory bandwidth; more
threads than pairs are never created.

Both screening choices use the same harmonized instruments:

- `IVW_FIXED` uses the fixed-effect IVW p-value;
- `IVW_MULTIPLICATIVE_RANDOM` uses the multiplicative-random IVW p-value and is
  the default in `XwasMrOptions.defaults(filter)`.

Only a threshold-passing pair runs MR-Egger, weighted median, bootstrap,
leave-one-out, instrument-strength warnings, and the remaining combined
analysis. Robust methods such as RAPS, contamination mixture, and the
PRESSO-style diagnostic can then be added for the much smaller hit list.

## Inspect results and accounting

```java
System.out.printf(
    "pairs=%d, analyzable=%d, hits=%d, below=%d, insufficient=%d, failures=%d%n",
    result.totalPairs(),
    result.analyzablePairs(),
    result.hits().size(),
    result.belowThresholdPairs(),
    result.insufficientInstrumentPairs(),
    result.failures().size());

for (XwasMrHit hit : result.hits()) {
    MrEstimate estimate = hit.screeningEstimate();
    System.out.printf("%s\t%s\t%s\t%g\t%g\t%g\t%g%n",
        hit.exposureId(), hit.outcomeCategory(), hit.outcomeId(),
        estimate.estimate(), estimate.standardError(), estimate.pValue(),
        hit.negativeLog10PValue());

    MrAnalysisResult diagnostics = hit.analysis();
    diagnostics.warnings().forEach(System.out::println);
    hit.harmonizationExclusions().forEach(System.out::println);
}
```

Pairs with fewer than three harmonized instruments are counted separately and
do not run the combined diagnostics. A pair can pass IVW screening but fail a
later diagnostic, for example when MR-Egger has no exposure-effect variation;
such cases appear in `failures()` with exposure/outcome IDs, exception type,
and message rather than stopping the whole scan.

The hit list and failures are deterministic for a fixed input order, options,
and seed. They remain in the same order with one thread or many threads.

## Follow-up analyses for hits

Apply sensitivity analysis to each retained hit rather than to every null cell:

```java
for (XwasMrHit hit : result.hits()) {
    List<HarmonizedInstrument> instruments = /* retain during input mapping */;
    MrRapsResult raps = RobustMendelianRandomization.raps(instruments);
    ContaminationMixtureResult mixture =
        ContaminationMixture.fit(instruments, 1001);
    MrPressoResult presso = MrPresso.analyze(instruments, 0.05);
    // Persist estimates, convergence, outliers, and diagnostic p-values.
}
```

`XwasMrHit.analysis()` already contains IVW, MR-Egger, weighted median, Wald
ratios, F statistics, heterogeneity, Egger intercept, leave-one-out, and
warnings. The pipeline retains harmonization exclusions but does not currently
retain the full harmonized instrument list in each hit to keep result memory
smaller. Applications needing additional robust calls should cache or rebuild
the small harmonized list for hits only.

Generate forest, funnel, scatter, and leave-one-out figures from the hit result
objects using the export schema in the
[end-to-end MR plotting section](mr-end-to-end.md#10-plotting).

## Performance and provenance checklist

- Clump once per exposure and record the reference manifest, population,
  window, r-squared threshold, and absent variants.
- Load and index each outcome once; do not scan a genome-wide outcome file for
  every molecular feature.
- Use bounded blocks rather than one future per pair.
- Screen with IVW before running bootstraps and leave-one-out diagnostics.
- Keep a fixed exposure order, outcome order, thread count, block size, and
  bootstrap seed in run metadata.
- Record total, analyzable, below-threshold, insufficient-instrument, failed,
  and retained pair counts.
- Correct across the full prespecified exposure-by-outcome testing family.
- Stratify summaries by the stored phenotype category without changing the
  original multiple-testing family after observing results.
- Reproduce the [xWAS MR benchmark](../xwas-mr-cli-performance.md) and profile
  representative molecular features before launching a much larger or
  materially different scan.

This design avoids redundant clumping and full diagnostics while preserving
the complete scientific two-level loop.
