# Mendelian randomization

> **v0.1.0 performance status:** This implementation is tested for numerical
> behavior but is not yet optimized or benchmarked as a high-throughput MR
> pipeline.

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
