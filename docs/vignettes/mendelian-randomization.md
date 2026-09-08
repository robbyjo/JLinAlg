# Mendelian randomization

> **v0.2.0 performance status:** Core estimators have independent numerical
> fixtures. The bounded xWAS screening path has a reproducible 45,000-pair
> Java/base-R benchmark; sensitivity paths remain workload-dependent.

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

Check `raps.converged()` before interpreting its estimate. This implementation
uses a Huber adjusted profile score, a nonnegative residual-moment estimate
of overdispersion, and a beta sandwich conditional on that plug-in dispersion.
The variance moment is truncated **after** averaging, not per instrument, and
the sandwich meat uses the same adjusted score as the estimating equation.
It is not the full jointly calibrated `mr.raps` procedure: the variance moment
is not outlier-robust and dispersion-estimation uncertainty is not included.

The contamination mixture is a bounded grid sensitivity model with valid
probabilities 0.05–0.95 and fixed additional invalid-effect variance 0.01 on
the outcome-effect scale. It is not R ConMix/MRCML. Its log-likelihood is
evaluated without probability floors; unresolved/boundary curvature throws
instead of substituting the grid spacing as an SE. Grid resolution and this
scale-dependent variance remain substantive modeling choices.

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

`direct.marginalFStatistics()` reports mean squared marginal exposure z
scores, not conditional instrument strength. The deprecated
`conditionalFStatistics()` now throws explicitly: the former values were
misnamed, and the input API has no specified cross-exposure sampling covariance
model for conditional F inference. Multivariable Egger uses the supplied
joint allele orientation; choose and record a common orientation convention
before calling it. Unlike univariate Egger, it does not automatically choose
an exposure-increasing reference trait.

Winner's-curse correction solves the exact two-sided truncated-normal
likelihood score. It rejects observations that fail the stated selection
event and thresholds above 1e6. Tail probabilities are evaluated in log space
and inverse Mills ratios by a stable tail continued fraction. For example,
an observed z=40.01 selected at 40 has corrected standardized effect
0.1071804266, not the old probability-floor/grid artifact 2.96074. This is an
effect correction, not an additional conditional SNP-association p-value.

Core IVW scales weights and variables before accumulation and reconstructs
the estimate and SE with separated binary exponents, so an overflowing
intermediate scale ratio does not erase an otherwise finite result. Univariate Egger
uses anchored weighted centering rather than subtracting nearly equal raw
moments. The regression suite includes an exact slope-2 line near exposure
offset 1e8 (formerly returned slope 4), IVW/Steiger unit scalings from 1e-200
to 1e200, allele-orientation tests, and the previously validated conditional,
generalized, and overlap-aware paths. These numerical fixes do not remove
the NOME, InSIDE, independence, or selection-model assumptions.

GRMs accept valid zero-variance/mean-imputed sample rows and singular PSD
matrices; supplied matrices are checked rather than repaired. PLINK clumping's
unphased haplotype likelihood is maximized over all admissible stationary
phases and both endpoints. The self-consistency cubic is partitioned at its
derivative roots, then each monotone interval is bisected; the observed
genotype likelihood selects the global candidate. There are no EM starts to
time out or stationary midpoints silently mistaken for maxima. Equivalent
phase maxima need not identify a unique haplotype phase even when r² is stable.

See [audit commands, reference checks, and raw warm timings](../../src/benchmark/resources/genetic-audit/AUDIT.md).

For estimator assumptions and limitations, see the
[MR scope document](../mr-timeseries-susie-sem.md) and
[numerical contract](../numerical-contract.md).
