# Mendelian randomization

> **v0.1.0 performance status:** This implementation is tested for numerical
> behavior but is not yet optimized or benchmarked as a high-throughput MR
> pipeline.

For a single workflow connecting CLI preparation through diagnostics and
figures, start with the [end-to-end MR vignette](mr-end-to-end.md).

## Find and prepare instruments from the CLI

Search the public NHGRI-EBI GWAS Catalog by trait. The search returns only
studies with downloadable full summary statistics and includes accession,
ontology trait, ancestry, sample description, license, and source URL:

```console
java -jar jlinalg-<version>.jar mr-instruments search --trait BMI --limit 20
```

`BMI` and a small set of common genetic-epidemiology acronyms are expanded to
their trait names before the ontology-backed search. Select a study after
checking phenotype definition, ancestry, sample composition, genome build,
and license, then stream its significant candidate variants to a local table:

```console
java -jar jlinalg-<version>.jar mr-instruments download \
  --study GCST... --out bmi-instruments.tsv --p-threshold 5e-8
```

For a local GWAS or molecular-QTL result, common column aliases are detected
automatically. Explicit mappings use canonical-target-to-source order:

```console
java -jar jlinalg-<version>.jar mr-instruments format \
  --input exposure.csv.gz --out exposure.mr.tsv --trait BMI \
  --map SNP=rsid,beta=estimate,se=stderr,eaf=frequency,\
    effect_allele=allele1,other_allele=allele2,pval=p_value
```

Output follows the MRInstruments/TwoSampleMR naming convention: `Phenotype`,
`SNP`, `beta`, `se`, `eaf`, `effect_allele`, `other_allele`, `pval`, `units`,
`ncase`, `ncontrol`, `samplesize`, and `gene`. Rows above the optional p-value
threshold or incompatible with the current biallelic-SNP MR layer are counted.
P-value filtering creates candidates, not independent instruments; perform
ancestry-matched LD clumping before using the independent-instrument methods:

```console
java -jar jlinalg-<version>.jar ld-db download \
  --database 1000g-phase3 --location /data/ld/1000g-phase3
java -jar jlinalg-<version>.jar clump \
  --database /data/ld/1000g-phase3 --population EUR \
  --instrument bmi-instruments.tsv --ld-threshold 0.001 \
  --output bmi-instruments.clumped.tsv
```

The default population (`EUR`), window (10,000 kb), r-squared threshold
(`0.001`), and index p-value threshold (`1`) match
`TwoSampleMR::clump_data`. Select a panel matching the exposure GWAS ancestry.
The command uses `pval.exposure`, `pval.outcome`, or `pval` automatically and
clumps each `id.exposure`/`id` group separately. It preserves every input column
for retained variants and writes retained rows in their original order.

## Harmonize exposure and outcome associations

Each summary association names the reported effect allele and other allele.
Effect-allele frequencies help resolve palindromic variants:

```java
List<SummaryAssociation> exposure = List.of(
    new SummaryAssociation("rs1", "A", "C", 0.10, 0.02, 0.20),
    new SummaryAssociation("rs2", "G", "T", 0.08, 0.02, 0.35),
    new SummaryAssociation("rs3", "C", "A", 0.12, 0.03, 0.40));

List<SummaryAssociation> outcome = List.of(
    new SummaryAssociation("rs1", "A", "C", 0.05, 0.02, 0.21),
    new SummaryAssociation("rs2", "T", "G", -0.03, 0.02, 0.65),
    new SummaryAssociation("rs3", "A", "C", -0.06, 0.03, 0.60));

HarmonizationResult aligned = AlleleHarmonizer.harmonize(
    exposure, outcome, HarmonizationOptions.defaults());
System.out.println(aligned.exclusions());
```

Never discard exclusions silently. Frequency-inconsistent, ambiguous
palindromic, duplicate, missing, and incompatible variants carry structured
reasons that belong in analysis provenance.

## Core independent-instrument analysis

At least three retained instruments are required because the combined result
includes MR-Egger and leave-one-out analyses:

```java
MrAnalysisResult mr = MendelianRandomization.analyze(
    aligned.instruments(), new MrOptions(0.95, 10_000, 20260901L));

MrEstimate fixedIvw = mr.ivwFixed();
MrEstimate randomIvw = mr.ivwMultiplicativeRandom();
MrEggerResult egger = mr.egger();

System.out.printf("IVW beta=%g se=%g p=%g Q=%g I2=%g%n",
    randomIvw.estimate(), randomIvw.standardError(),
    randomIvw.pValue(), randomIvw.cochranQ(),
    randomIvw.iSquared());
```

Also inspect the MR-Egger intercept, I-squared GX, mean instrument F statistic,
leave-one-out estimates, and warnings. Sensitivity estimators do not turn the
exclusion restriction or InSIDE assumption into testable facts.

## LD-aware generalized IVW and Egger

Supply an allele-aligned positive-definite correlation matrix in exactly the
harmonized variant order:

```java
double[][] ld = {
    {1.00, 0.15, 0.02},
    {0.15, 1.00, 0.10},
    {0.02, 0.10, 1.00}
};

CorrelatedMrEstimate generalized =
    CorrelatedMendelianRandomization.ivw(
        aligned.instruments(), ld, true);
CorrelatedMrEggerResult generalizedEgger =
    CorrelatedMendelianRandomization.egger(
        aligned.instruments(), ld);
```

The sign of each LD row/column must correspond to that instrument's aligned
exposure effect allele. Prune or scientifically regularize a nearly singular
matrix before calling this API.

## Directionality and robustness

```java
SteigerResult steiger = SteigerFiltering.analyze(
    aligned.instruments(), exposureSampleSizes, outcomeSampleSizes);

MrRapsResult raps = RobustMendelianRandomization.raps(
    aligned.instruments());
MrPressoResult presso = MrPresso.analyze(
    aligned.instruments(), 0.05);
ContaminationMixtureResult mixture = ContaminationMixture.fit(
    aligned.instruments(), 1001);
```

The PRESSO-style result is an analytic robust-center/outlier diagnostic rather
than the simulation calibration of the R MR-PRESSO package. Report which
variants were removed and show estimates before and after exclusion.

## Multivariable and overlapping-sample MR

```java
MultivariableMrResult direct =
    MultivariableMendelianRandomization.ivw(
        multivariableInstruments,
        List.of("LDL", "HDL"));

OverlapAwareMrResult overlapAware =
    OverlapAwareMendelianRandomization.ivw(
        aligned.instruments(), exposureOutcomeCovariance);
```

`exposureOutcomeCovariance` is per-instrument sampling covariance, not a sample
overlap percentage. For selection sensitivity, apply
`WinnerCurseCorrection.correct(effect, se, selectionZ)` to the appropriate
exposure associations and retain both corrected and original results.

For estimator assumptions and limitations, see the
[MR scope document](../mr-timeseries-susie-sem.md) and
[numerical contract](../numerical-contract.md).
