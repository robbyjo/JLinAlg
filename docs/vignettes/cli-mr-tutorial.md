# CLI-only Mendelian randomization: download to results

This page is a file-by-file workflow for users who do not write Java. It
covers public GWAS discovery, custom instrument files, normalization, LD
clumping, allele harmonization inside JLinAlg, causal estimation, and output
files. Replace `jlinalg-<version>.jar` with the downloaded executable JAR.

## 1. Decide the analysis before downloading data

Write down the exposure, outcome, units, ancestry, genome build, covariate
adjustments, sample sizes, and likely sample overlap. Exposure p-values select
candidate instruments; outcome p-values must not be used to choose them.
Prefer an LD reference matched to the exposure GWAS ancestry. Confirm that
the summary-statistic license permits the intended use.

## 2. Find a suitable public exposure GWAS

JLinAlg searches the public NHGRI-EBI GWAS Catalog:

```powershell
java -jar jlinalg-<version>.jar mr-instruments search --trait "body mass index" --limit 20
```

Review the accession, exact phenotype definition, ancestry, sample
description, genome build, and availability of full summary statistics.
Choose an accession deliberately; the first search hit is not automatically
the best scientific match.

Download genome-wide-significant candidates:

```powershell
java -jar jlinalg-<version>.jar mr-instruments download `
  --study GCST00000000 --out bmi-candidates.tsv --p-threshold 5e-8
```

Replace the fictitious accession. Downloading selects by p-value but does not
make variants independent; LD clumping is a separate step.

## 3. Use your own candidate-instrument file

A usable instrument row needs the variant ID, exposure effect, standard
error, effect allele, other allele, and p-value. Effect-allele frequency is
optional but strongly recommended for resolving palindromic alleles. A trait
label is required when several exposures share a file.

JLinAlg's normalized columns are:

```text
Phenotype SNP beta se eaf effect_allele other_allele pval
units ncase ncontrol samplesize gene
```

CSV, TSV, `.csv.gz`, and `.tsv.gz` are accepted. Common names such as
`rsid`, `effect_size`, `stderr`, `a1`, `a2`, and `p_value` are
detected. For a study-specific table, map canonical target names to source
names:

```powershell
java -jar jlinalg-<version>.jar mr-instruments format `
  --input my-exposure.csv.gz --out exposure-candidates.tsv `
  --trait "LDL cholesterol" --p-threshold 5e-8 `
  --map SNP=MarkerName,beta=Effect,se=StdErr,eaf=Freq1,effect_allele=Allele1,other_allele=Allele2,pval=Pvalue
```

The syntax is `TARGET=SOURCE`: do not rename the source file by hand unless
you prefer to. `--trait` supplies `Phenotype` when the input lacks a trait
column. If Effect contains odds ratios or hazard ratios rather than additive
effects, add `--effect-scale odds-ratio` or
`--effect-scale hazard-ratio`; JLinAlg logs the ratios for MR.

A file containing only chosen rsIDs is not sufficient for MR. Join those IDs
to their exposure beta, SE, alleles, and preferably EAF first, then run the
formatter. Invalid or non-biallelic rows are counted and omitted.

## 4. Install and select an LD reference

```powershell
java -jar jlinalg-<version>.jar ld-db list
java -jar jlinalg-<version>.jar ld-db download `
  --database 1000g-phase3 --location D:\reference\1000g-phase3
```

The 1000 Genomes Phase 3 package provides AFR, AMR, EAS, EUR, and SAS panels.
Select the panel matching the exposure study rather than accepting EUR merely
because it is the default. Keep the database manifest with analysis
provenance.

## 5. LD-clump the exposure

```powershell
java -jar jlinalg-<version>.jar clump `
  --database D:\reference\1000g-phase3 --population EUR `
  --instrument exposure-candidates.tsv `
  --clump-kb 10000 --ld-threshold 0.001 --p-threshold 1 `
  --output exposure-clumped.tsv
```

`--p-threshold 1` tells clumping to use every already-selected candidate.
The defaults are a 10,000 kb window, r-squared 0.001, and EUR. The output
retains the original columns and input order. Variants absent from the
reference are excluded; a large absent fraction suggests an identifier,
genome-build, or ancestry mismatch.

The clumper recognizes `SNP` or `rsid` and common p-value names. Override
unusual files with `--snp-column`, `--pval-column`, and `--id-column`.
When several exposures occupy one table, each exposure ID group is clumped
independently.

## 6. Prepare the outcome summary statistics

The outcome must contain associations for the clumped exposure variants, with
beta, SE, effect allele, other allele, and preferably EAF. Normalize a local
outcome table without an outcome-significance filter:

```powershell
java -jar jlinalg-<version>.jar mr-instruments format `
  --input coronary-disease.tsv.gz --out outcome-mr.tsv `
  --trait "coronary artery disease" --p-threshold 1 `
  --map SNP=variant,beta=log_odds,se=standard_error,eaf=effect_af,effect_allele=effect_allele,other_allele=other_allele,pval=p_value
```

If an outcome service can query only requested variants, query the SNP column
of `exposure-clumped.tsv`. Resolve genome builds before MR; allele
harmonization is not liftover.

## 7. Run one exposure-outcome MR entirely from the CLI

The `mr-xwas` command also handles a single exposure and single outcome. It
groups the files, intersects variants, harmonizes alleles, and estimates MR.
Use `--p-threshold 1` to retain the one analyzable pair regardless of its
causal p-value:

```powershell
java -jar jlinalg-<version>.jar mr-xwas `
  --exposure exposure-clumped.tsv --outcome outcome-mr.tsv `
  --output bmi-cad-mr.tsv --p-threshold 1 `
  --screen-method random --threads 1 `
  --fdr-output bmi-cad-all-pairs.tsv `
  --follow-up-output bmi-cad-follow-up.tsv
```

Required exposure columns are `SNP`, `beta`, `se`,
`effect_allele`, `other_allele`, and one of `id.exposure`, `gene`,
or `Phenotype`. Required outcome columns are the equivalent fields and
`id.outcome` or `Phenotype`. Use `--exposure-id-column` and
`--outcome-id-column` for custom group columns. EAF is optional. Palindromic
variants without enough frequency information may be excluded rather than
guessed.

For many exposures and outcomes, concatenate normalized rows into long-format
tables and use more threads. Clump each exposure once, not once per outcome.
Set exactly one screening threshold: `--p-threshold`,
`--log10-p-threshold`, or `--negative-log10-p-threshold`.

## 8. Know every output file

For the command above:

- `bmi-cad-mr.tsv` contains retained pairs and the main IVW estimate, SE,
  causal p-value, confidence interval, heterogeneity, mean F, MR-Egger,
  weighted median, harmonization exclusions, and warnings.
- `bmi-cad-mr.failures.tsv` is always written. It records pair identifiers,
  exception type, and message for diagnostic failures.
- `bmi-cad-all-pairs.tsv` contains every successfully screened pair plus
  `threshold_passed` and BH `fdr_bh`. BH uses the complete analyzable
  family, not only retained hits.
- `bmi-cad-follow-up.tsv` is optional and contains hit-only RAPS,
  contamination-mixture, and PRESSO-style results.

Pairs with fewer than three harmonized instruments have no full diagnostic
fit. Inspect exclusion counts and the failures file before treating an empty
result as a null finding.

## 9. Analyze an already harmonized table

If another tool has already aligned outcome effects to the exposure effect
allele, prepare:

```text
variant_id beta_exposure se_exposure beta_outcome se_outcome
rs101      0.08          0.02        0.03         0.01
```

Then run:

```powershell
java -jar jlinalg-<version>.jar mr-estimate --input harmonized.tsv `
  --method all --output mr-estimates.tsv --plot mr-scatter.svg
```

The estimator output columns are `method`, `estimate`,
`standard_error`, `causal_p_value`, `ci_lower`, `ci_upper`, and
`instruments`. Methods include `ivw-fixed`, `ivw-random`, `egger`,
`weighted-median`, generalized LD-aware variants, overlap-aware IVW,
conditional association, and `all`. A signed, allele-aligned LD matrix
passed with `--ld` must follow the exact input row order.

`mr-estimate` assumes harmonization; it has no allele columns and cannot
detect an allele flip. If you have separate unharmonized exposure and outcome
files, use `mr-xwas` instead.

## 10. Interpretation checklist

- Treat strength, heterogeneity, Egger intercept, robust estimates, and
  leave-one-out behavior as a joint diagnostic set.
- Record all harmonization exclusions and variants missing from the LD panel.
- Exponentiate an estimate only when its outcome scale is log odds or log
  hazard and state the exposure unit.
- Reverse-direction MR is a second analysis with independently selected
  instruments, not a relabeling of the first result.
- Colocalization tests whether regional signals appear shared; it does not
  turn an MR association into proof of causality.

Continue with [CLI-only SuSiE](cli-susie-tutorial.md) and
[CLI-only colocalization](cli-colocalization-tutorial.md) for locus follow-up.
The [MR estimator vignette](mendelian-randomization.md) provides the statistical
details behind these commands.
