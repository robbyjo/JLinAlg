# End-to-end Mendelian randomization

This vignette joins JLinAlg's command-line data preparation to its Java MR
API. It covers instrument discovery or custom QTL input, ancestry-matched LD
clumping, harmonization, estimation, diagnostics, plotting, bidirectional MR,
and molecular-trait studies.

> **Current interface boundary:** `ld-db`, `mr-instruments`, and `clump` are
> executable CLI commands. MR estimation is currently a Java API; there is no
> `mr analyze` file-to-file command yet. JLinAlg also exposes plot-ready result
> objects but does not currently render plots. The plotting section below gives
> an explicit export schema and R examples so these boundaries are reproducible
> rather than hidden.

## Workflow at a glance

```text
exposure GWAS/QTL -> format/select candidates -> ancestry-matched clumping
                                                   |
outcome GWAS -----------------------> allele harmonization
                                                   |
                            IVW/Egger/median + diagnostics
                                                   |
                           tables, plots, sensitivity analyses
```

Choose the exposure and outcome datasets before running commands. Record their
phenotype definitions, units, genome builds, ancestries, sample sizes, sample
overlap, covariate adjustments, and licenses. The exposure associations select
instruments; the outcome associations must provide results for those variants.
Do not select instruments using the outcome p-values.

## 1. Install an LD reference

List the available databases and install one into a durable directory:

```console
java -jar jlinalg-<version>.jar ld-db list
java -jar jlinalg-<version>.jar ld-db download \
  --database 1000g-phase3 --location /data/ld/1000g-phase3
```

The installation contains AFR, AMR, EAS, EUR, and SAS PLINK panels when the
source provides them. Select the panel that best matches the exposure GWAS.
The default is `EUR`; a default is not evidence that EUR is appropriate.

## 2A. Find public exposure instruments

Search the NHGRI-EBI GWAS Catalog by trait. Common acronyms such as `BMI` are
expanded before the ontology-backed search:

```console
java -jar jlinalg-<version>.jar mr-instruments search \
  --trait BMI --limit 20
```

Inspect the returned accession, exact trait, ancestry, sample description,
genome build, and license. Then download significant candidate variants:

```console
java -jar jlinalg-<version>.jar mr-instruments download \
  --study GCST... --out bmi.candidates.tsv --p-threshold 5e-8
```

This p-value filter creates candidate instruments. It does not establish
independence, strength in the analysis sample, or the exclusion restriction.

## 2B. Format custom GWAS or QTL instruments

For a local GWAS, eQTL, mQTL, pQTL, or metabolite-QTL table, normalize common
columns automatically or map study-specific names explicitly:

```console
java -jar jlinalg-<version>.jar mr-instruments format \
  --input exposure.csv.gz --out exposure.candidates.tsv \
  --trait BMI --p-threshold 5e-8 \
  --map SNP=rsid,beta=estimate,se=stderr,eaf=frequency,\
effect_allele=allele1,other_allele=allele2,pval=p_value
```

Canonical output columns are `Phenotype`, `SNP`, `beta`, `se`, `eaf`,
`effect_allele`, `other_allele`, `pval`, `units`, `ncase`, `ncontrol`,
`samplesize`, and `gene`. Use `--effect-scale odds-ratio` or
`--effect-scale hazard-ratio` when the mapped effect is a ratio; JLinAlg stores
the corresponding log effect for MR.

When a file contains several molecular features, supply an `id.exposure`
column before clumping, or map an equivalent column with `--id-column` during
clumping. Every feature is then clumped independently. Without an ID column,
all rows are treated as one exposure.

## 3. LD-clump the exposure candidates

```console
java -jar jlinalg-<version>.jar clump \
  --database /data/ld/1000g-phase3 \
  --population EUR \
  --instrument bmi.candidates.tsv \
  --ld-threshold 0.001 \
  --output bmi.instruments.tsv
```

The defaults match `TwoSampleMR::clump_data`: a 10,000 kb window, r-squared
threshold 0.001, index p-value threshold 1, and EUR population. The command
uses `pval.exposure`, `pval.outcome`, `pval`, or `p` in that order and clumps
each `id.exposure` or `id` separately. Use `--clump-kb`, `--p-threshold`,
`--snp-column`, `--pval-column`, and `--id-column` when needed.

The output retains every original column and preserves input order. For a
multi-variant group, variants absent from the reference are excluded. Review
how many variants were absent; systematic loss can indicate an ancestry,
genome-build, or identifier mismatch.

## 4. Prepare outcome associations

Obtain outcome effects for every retained exposure instrument. A local outcome
table can use the same formatter, but do not apply an exposure-style
significance threshold:

```console
java -jar jlinalg-<version>.jar mr-instruments format \
  --input outcome.tsv.gz --out outcome.mr.tsv --trait CHD \
  --map SNP=variant,beta=log_odds,se=standard_error,eaf=eaf,\
effect_allele=effect_allele,other_allele=other_allele,pval=p_value
```

Subset or query the outcome by the SNPs in `bmi.instruments.tsv`. For
two-sample MR, prefer non-overlapping samples when practical. If samples
overlap, record or estimate the per-variant covariance between exposure and
outcome estimates; a sample-overlap percentage alone is not the covariance
required by JLinAlg's overlap-aware estimator.

## 5. Map rows into the Java API and harmonize alleles

For each normalized row, construct a `SummaryAssociation` as follows:

```java
SummaryAssociation association = new SummaryAssociation(
    row.get("SNP"),
    row.get("effect_allele"),
    row.get("other_allele"),
    Double.parseDouble(row.get("beta")),
    Double.parseDouble(row.get("se")),
    Double.parseDouble(row.get("eaf")));
```

Build `List<SummaryAssociation> exposure` from the clumped exposure file and
`List<SummaryAssociation> outcome` from the outcome file, then harmonize:

```java
HarmonizationResult harmonized = AlleleHarmonizer.harmonize(
    exposure, outcome, HarmonizationOptions.defaults());

List<HarmonizedInstrument> instruments = harmonized.instruments();
harmonized.exclusions().forEach(System.out::println);
```

Harmonization aligns outcome effects to the exposure effect allele, including
allele swaps and strand complements. Palindromic variants are retained only
when the two effect-allele frequencies select an unambiguous orientation.
Always retain the structured exclusions as analysis provenance.

JLinAlg currently accepts biallelic SNP alleles `A`, `C`, `G`, and `T`, finite
effects, and positive finite standard errors. Resolve genome builds and variant
identifiers before this step; harmonization does not perform liftover.

## 6. Run the primary MR analysis

The combined analysis requires at least three harmonized instruments because
it includes MR-Egger, weighted median, and leave-one-out results:

```java
MrOptions options = new MrOptions(0.95, 10_000, 20260904L);
MrAnalysisResult mr = MendelianRandomization.analyze(instruments, options);

MrEstimate fixed = mr.ivwFixed();
MrEstimate random = mr.ivwMultiplicativeRandom();
MrEggerResult egger = mr.egger();
MrEstimate median = mr.weightedMedian();

System.out.printf("IVW random: beta=%g, SE=%g, P=%g, 95%% CI=[%g, %g]%n",
    random.estimate(), random.standardError(), random.pValue(),
    random.confidenceLower(), random.confidenceUpper());
mr.warnings().forEach(System.out::println);
```

The primary result should normally be chosen from the scientific design, not
after comparing p-values. JLinAlg provides fixed IVW, multiplicative-random
IVW, MR-Egger, weighted median, and per-variant Wald ratios. Binary-outcome
effects are on the log-odds scale when that was the outcome GWAS scale;
`Math.exp(estimate)` converts a log-odds causal estimate to an odds ratio.

One instrument can be analyzed with
`MendelianRandomization.waldRatios(instruments)`. Two instruments can use
`MendelianRandomization.ivw(...)`. The combined `analyze` method intentionally
requires three.

## 7. Diagnostics and sensitivity analyses

No single diagnostic proves that all instruments are valid. Interpret the
following together with biological annotation and study design.

| Question | JLinAlg output or call | Interpretation |
| --- | --- | --- |
| Instrument strength | `waldRatios().fStatistic()` and `meanFStatistic()` | Approximate per-SNP and mean `(betaX/seX)^2`; values below 10 trigger a warning, but a threshold is not a guarantee against weak-instrument bias. |
| IVW heterogeneity | `ivwFixed().cochranQ()`, `heterogeneityDegreesOfFreedom()`, `heterogeneityPValue()`, `iSquared()` | Excess dispersion can reflect invalid instruments, different mechanisms, or model misspecification. Random IVW inflates uncertainty by `max(1, Q/df)`. |
| Egger heterogeneity | The same fields on `egger().slope()` | Residual heterogeneity around the Egger line. |
| Directional pleiotropy | `egger().intercept()` and `interceptPValue()` | An intercept different from zero is evidence compatible with directional horizontal pleiotropy under Egger assumptions; low power is common. |
| Egger reliability | `egger().iSquaredGx()` | Low I-squared GX warns about dilution from imprecise exposure effects. |
| Influential SNP | `leaveOneOut()` | Repeats fixed IVW after omitting each instrument; look for a materially shifted estimate or interval. |
| Directionality | `SteigerFiltering.analyze(...)` | Compares estimated exposure and outcome variance explained per SNP and in aggregate. It requires exposure and outcome sample size for every SNP. |
| Robustness/outliers | weighted median, RAPS, contamination mixture, and `MrPresso.analyze(...)` | Sensitivity to invalid or outlying instruments under different assumptions. |

Example diagnostic calls:

```java
System.out.printf("IVW Q=%g, df=%d, P=%g, I2=%g%n",
    fixed.cochranQ(), fixed.heterogeneityDegreesOfFreedom(),
    fixed.heterogeneityPValue(), fixed.iSquared());

System.out.printf("Egger intercept=%g, SE=%g, P=%g, I2GX=%g%n",
    egger.intercept(), egger.interceptStandardError(),
    egger.interceptPValue(), egger.iSquaredGx());

mr.leaveOneOut().forEach(value -> System.out.printf(
    "omit %s: beta=%g, SE=%g, P=%g%n",
    value.omittedVariantId(), value.estimate().estimate(),
    value.estimate().standardError(), value.estimate().pValue()));

SteigerResult steiger = SteigerFiltering.analyze(
    instruments, exposureSampleSizes, outcomeSampleSizes);

MrRapsResult raps = RobustMendelianRandomization.raps(instruments);
ContaminationMixtureResult mixture = ContaminationMixture.fit(instruments, 1001);
MrPressoResult presso = MrPresso.analyze(instruments, 0.05);
```

Check `raps.converged()`. JLinAlg's PRESSO-style procedure is a fast analytic
robust-center/outlier diagnostic; it is not the simulation-calibrated R
MR-PRESSO algorithm. Report the raw and corrected results, global p-value,
outliers, and distortion instead of silently deleting variants.

Steiger directionality is a sensitivity analysis, not a replacement for a
reverse-direction MR. Measurement error, trait scaling, ascertainment, and
sample-size differences can affect its comparison.

## 8. Bidirectional (two-way) MR

There is no special bidirectional wrapper because bidirectional MR consists of
two separately designed analyses:

1. Select and clump instruments for `X`; obtain their associations with `Y`;
   harmonize and run `X -> Y`.
2. Independently select and clump instruments for `Y`; obtain their
   associations with `X`; harmonize and run `Y -> X`.

Do not obtain the reverse result by merely swapping beta columns for the same
instrument set. The reverse analysis needs instruments selected for the new
exposure, with its own strength, LD, harmonization, and sensitivity checks.
Call `MendelianRandomization.analyze` once for each harmonized direction and
report both instrument sets. Bidirectional MR is also distinct from
multivariable MR.

## 9. Multivariable, correlated, and overlapping-sample MR

For multiple proposed exposures measured against the same outcome, construct
`MultivariableInstrument` rows and fit direct effects:

```java
MultivariableMrResult mvivw =
    MultivariableMendelianRandomization.ivw(values, exposureNames);
MultivariableMrResult mvegger =
    MultivariableMendelianRandomization.egger(values, exposureNames);
```

For intentionally retained correlated cis variants, supply an allele-aligned,
positive-definite LD correlation matrix in the exact instrument order:

```java
CorrelatedMrEstimate givw =
    CorrelatedMendelianRandomization.ivw(instruments, alleleAlignedLd, true);
CorrelatedMrEggerResult gegger =
    CorrelatedMendelianRandomization.egger(instruments, alleleAlignedLd);
```

The `clump` command produces independent index variants; it does not currently
export the signed LD matrix needed by these generalized estimators.

When exposure and outcome estimates have overlapping samples:

```java
OverlapAwareMrResult overlap =
    OverlapAwareMendelianRandomization.ivw(instruments, samplingCovariance);
```

`samplingCovariance[j]` is the covariance of the exposure and outcome effect
estimates for SNP `j`, not an overlap fraction. Check `overlap.converged()`.

## 10. Plotting

JLinAlg does not yet include a plotting dependency or SVG/PNG MR renderer. The
result objects expose everything needed for the conventional plots:

| Plot | Source data |
| --- | --- |
| SNP forest | `mr.waldRatios()`: variant, estimate, standard error; calculate `estimate +/- 1.959964 * SE` |
| Funnel | The same Wald ratios: x = estimate, y = `1/SE`; add the selected pooled estimate as a vertical line |
| Leave-one-out | `mr.leaveOneOut()`: omitted variant and nested estimate/CI |
| Exposure-outcome scatter | `HarmonizedInstrument`: exposure effect/SE and outcome effect/SE; overlay IVW and Egger lines |
| Method forest | `ivwFixed`, `ivwMultiplicativeRandom`, `egger().slope()`, `weightedMedian`, and optional robust estimates |

For example, export `wald.tsv` with columns `variant`, `estimate`, `se`,
`lower`, and `upper`, and `leave_one_out.tsv` with `omitted`, `estimate`,
`lower`, and `upper`. These R commands render standard plots:

```r
library(ggplot2)

wald <- read.delim("wald.tsv")
ggplot(wald, aes(y = reorder(variant, estimate), x = estimate)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y") +
  geom_point() + labs(x = "Wald ratio (95% CI)", y = NULL)

ggplot(wald, aes(x = estimate, y = 1 / se)) +
  geom_vline(xintercept = ivw_estimate, linetype = 2) +
  geom_point() + labs(x = "Wald ratio", y = "Precision (1/SE)")

loo <- read.delim("leave_one_out.tsv")
ggplot(loo, aes(y = reorder(omitted, estimate), x = estimate)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y") +
  geom_point() + labs(x = "IVW estimate after omission", y = NULL)
```

Set `ivw_estimate` to the selected IVW estimate before drawing the funnel plot.
A funnel plot is a visual asymmetry diagnostic, not a formal pleiotropy test.
Preserve the numeric diagnostic tables alongside every figure.

## 11. TWAS, EWAS, PWAS, and metabolome-wide settings

The same summary-data MR estimators can be used when the exposure is a
molecular phenotype, but the input must still be variant-exposure QTL effects:

| Analysis label | Typical exposure association |
| --- | --- |
| transcript/proposed TWAS follow-up | eQTL effect on expression or splicing |
| EWAS follow-up | mQTL effect on DNA methylation at a CpG or region |
| PWAS/proteome-wide MR | pQTL effect on protein abundance |
| metabolome-wide MR | metabolite-QTL effect on metabolite abundance |

A TWAS gene-level association, EWAS CpG-outcome association, predicted protein
score, or metabolome-wide feature p-value is not by itself a genetic
instrument. Retain the underlying SNP-feature beta, SE, effect allele, other
allele, frequency, sample size, tissue/assay, and units.

Recommended design points for high-throughput molecular MR:

- Put the gene, transcript, CpG, protein, or metabolite identifier in
  `id.exposure` so clumping occurs separately per feature.
- Decide and report whether instruments are cis-only or include trans-QTLs.
  Trans instruments may be stronger but often create additional pleiotropy
  concerns.
- Match tissue, ancestry, genome build, allele coding, assay, and molecular
  trait transformation. A causal estimate is per unit of the reported QTL
  scale.
- Account for testing many features with a declared FDR or family-wise error
  procedure; do not reinterpret nominal MR p-values as study-wide evidence.
- For cis molecular MR, investigate colocalization or fine-mapping. LD clumping
  alone does not establish that the molecular trait and outcome share a causal
  variant. JLinAlg does not currently implement a colocalization model.
- Consider correlated cis-MR rather than discarding all but one signal when a
  scientifically justified allele-aligned LD matrix is available.
- Record sample overlap, winner's-curse risk, weak instruments, and tissue
  specificity. `WinnerCurseCorrection.correct(...)` is available for a
  selected normal association, but selection-aware study design remains
  preferable.
- Use multivariable MR when the scientific question concerns direct effects of
  several molecular traits; it is not a substitute for colocalization or
  mediation analysis.

For thousands of features, stream and partition the formatted tables by
`id.exposure`, keep a machine-readable exclusions table, and apply the same
predeclared pipeline to every feature. The v0.2.0 bounded xWAS screen is
numerically checked against base R and benchmarked for a high-throughput
molecular scan; hit-specific diagnostics remain workload-dependent.

## Reporting checklist

At minimum, preserve and report:

- exposure and outcome dataset accessions, versions, ancestry, units, genome
  build, sample sizes, covariates, and sample overlap;
- candidate selection threshold, LD database/version/population, clumping
  window and r-squared threshold, absent variants, and retained SNPs;
- harmonization exclusions and treatment of palindromic variants;
- per-SNP effects, SEs, alleles, frequencies, F statistics, and Wald ratios;
- prespecified primary estimator and every sensitivity estimator attempted;
- Q, degrees of freedom, heterogeneity p-value, I-squared, Egger intercept,
  I-squared GX, Steiger result, leave-one-out results, and robust/outlier checks;
- bidirectional instrument sets, if used, and multiple-testing correction for
  high-throughput molecular analyses;
- code version, random seed, database manifest, and warnings.

The [TwoSampleMR analysis vignette](https://github.com/MRCIEU/TwoSampleMR/blob/master/vignettes/perform_mr.Rmd)
provides a useful comparison for standard result tables and plots. Study
planning and reporting should also follow the
[2023 MR investigation guidelines](https://pmc.ncbi.nlm.nih.gov/articles/PMC7384151/)
and the [STROBE-MR explanation and checklist](https://pmc.ncbi.nlm.nih.gov/articles/PMC8546498/).
For molecular cis-MR, the guidelines emphasize that colocalization provides
complementary evidence and should be considered explicitly.

For estimator definitions and numerical assumptions, see the
[MR numerical contract](../numerical-contract.md) and
[MR scope document](../mr-timeseries-susie-sem.md).
