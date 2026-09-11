# CLI association tutorial: from one model to xWAS

This tutorial assumes only Java and the executable JLinAlg JAR. It begins with
one phenotype model, then adds numeric omics, genotype files, relatedness,
non-Gaussian outcomes, survival analysis, and penalized regression. Replace
`jlinalg-<version>.jar` with the downloaded filename.

Start any new analysis with `--dry-run`. It validates files and the formula
and prints the chosen model without fitting. Add `--explain` to print the
same routing information and continue with the analysis.

## 1. Begin with a phenotype file

A phenotype file is CSV or TSV with one observation per row. The column passed
to `--id` must be unique because it identifies that row and, when omics are
added later, matches the sample headers in the omics file.

```text
sample_id,BMI,Sex,Age,clinic
S001,24.8,F,52,A
S002,31.1,M,60,A
S003,27.5,F,47,B
S004,29.0,M,55,B
```

### Phenotype-only OLS

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --id sample_id `
  --formula "BMI ~ Sex + Age" --model ols --out bmi-ols.csv
```

The result has one row for every fitted coefficient. A `.csv` output is
comma-separated; `.tsv` and other suffixes are tab-separated. Rows missing
BMI, Sex, or Age are omitted as complete cases and counted on screen and in
`bmi-ols.csv.log`. The JSON manifest is `bmi-ols.csv.manifest.json`.

### Phenotype-only REML

Add a random term when observations are clustered. Here, clinic is an
ordinary independent random intercept:

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --id sample_id `
  --formula "BMI ~ Sex + Age + (1|clinic)" --model lmm `
  --variance-components refit --df satterth --out bmi-clinic-reml.csv
```

This is a likelihood-based Gaussian mixed model with variance components
estimated by REML. It is not P3D or EMMAX. `--df kr` requests
Kenward-Roger inference; `--df satterth` requests Satterthwaite inference.

## 2. Add a numeric omics matrix

Expression, methylation, protein, and generic numeric inputs are
feature-by-sample. The first column is the feature ID; every remaining column
header is a sample ID.

```text
gene_id,S001,S002,S003,S004
GENE_A,5.2,5.7,4.9,5.1
GENE_B,2.0,2.4,2.2,NA
```

JLinAlg intersects these headers with the phenotype `--id` values and
reorders automatically. It reports omics, phenotype, aligned, omitted, and
analysis sample counts. Numeric omics may be CSV or TSV.

### Omics OLS

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics expression.csv `
  --id sample_id --omics-type expression `
  --formula "BMI ~ Sex + Age + <omics>" --model ols `
  --transform "<omics> = winsor_mad(k=4) | zscore()" `
  --threads 8 --block-size auto --out bmi-expression-ols.csv
```

Each feature replaces `<omics>` in turn. Filtering and sample alignment occur
before row-wise transformation; missing feature values are mean-imputed for
the fit. Use `--annot genes.tsv --annot-id gene_id --annot-cols symbol,chr`
to append selected annotation columns.

### Omics REML

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics expression.csv `
  --id sample_id --omics-type expression `
  --formula "BMI ~ Sex + Age + <omics> + (1|clinic)" `
  --model lmm --variance-components refit --df satterth `
  --transform "<omics> = winsor_mad(k=4)" `
  --threads 8 --block-size auto --out bmi-expression-reml.csv
```

For non-genotype omics, `--variance-components auto` also chooses an exact
per-feature REML refit. The random-effects pathway is lmer-like, not a
single-null-model approximation.

## 3. Common-variant GWAS

Accepted genotype inputs are VCF, BGZF-compressed VCF, BCF, BGEN layout 2, and
additive-dosage CSV/TSV. VCF/BCF uses FORMAT/DS when available and otherwise
called GT. Multiallelic VCF records become one row per alternate allele. BGEN
supports zlib, zstd, or uncompressed layout-2 blocks; use
`--sample-file cohort.sample` when sample IDs are external.

A delimited dosage file can be compact:

```text
id,S001,S002,S003,S004
rs100,0,1,2,1
```

or carry genomic fields:

```text
id,chromosome,position,ref,alt,S001,S002,S003,S004
rs100,1,101001,A,G,0,1,2,1
```

Dosages are alternate-allele counts in [0,2]. Blank, `.`, `NA`, `N/A`,
`null`, and `NaN` are missing. For an ambiguously named delimited file,
specify `--omics-type gwas`.

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics cohort.vcf.gz `
  --id sample_id --formula "BMI ~ Sex + Age + PC1 + PC2 + <omics>" `
  --min-maf 0.01 --min-mac 20 --max-marker-missing 0.02 --min-info 0.8 `
  --threads 16 --block-size auto --out bmi-gwas.csv
```

`--min-info` requires an available imputation-quality value. Every excluded
variant remains in the audit table with `status=filtered` and a
`filter_reason`; only successful tests enter BH FDR. Add a labeled dense GRM
or GCTA prefix with `--grm cohort --individual-id subject_id` for a
genotype LMM. That high-throughput genotype path fits one REML null model and
uses P3D/EMMAX score tests, unlike numeric-omics REML refits.

## 4. Rare-variant GWAS

For a single-marker rare-variant scan, put an upper bound on frequency or
count. This example keeps variants with MAF at most 1% and MAC from 3 through
40:

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics cohort.bgen `
  --sample-file cohort.sample --id sample_id --omics-type gwas `
  --formula "BMI ~ Sex + Age + PC1 + PC2 + <omics>" `
  --max-maf 0.01 --min-mac 3 --max-mac 40 `
  --max-marker-missing 0.02 --out bmi-rare-single-variant.tsv
```

Frequency and count thresholds are applied to the aligned analysis cohort.
Single-marker tests can be poorly powered for very rare alleles. JLinAlg also
implements Burden, SKAT, and SKAT-O set tests in the Java API, but a
gene/set-membership CLI is not yet exposed; do not describe the command above
as SKAT or a gene-based test.

## 5. Pedigree files, slowly and explicitly

There can be three different identifiers:

1. The **omics ID** labels a measured sample, assay, or observation.
2. The **subject ID** labels a person and may repeat across visits or assays.
3. The **pedigree ID** labels that subject in the ancestry file.

Suppose the phenotype file contains both omics and subject IDs:

```text
omics_id,subject_id,BMI,Sex,Age,visit
OM9001,1001,24.8,F,52,1
OM9002,1002,31.1,M,60,1
OM9003,1003,27.5,F,28,1
OM9004,1003,28.0,F,30,2
OM9005,9009,26.2,M,44,1
```

The omics matrix headers are `OM9001` through `OM9005`, so use
`--id omics_id`. The pedigree uses subject IDs:

```text
family_id,member_id,parent1_id,parent2_id
F10,1001,0,0
F10,1002,0,0
F10,1003,1001,1002
```

Tell JLinAlg what each pedigree column means:

- `--pedigree-id member_id`: the individual represented by this row.
- `--sire-id parent1_id`: parent 1; the historical option name is retained.
- `--dam-id parent2_id`: parent 2.
- `--pedigree-family-id family_id`: optional family qualifier.
- `--individual-id subject_id`: phenotype column matched to pedigree members.

The formula pedigree term must use the phenotype matching-column name,
`(1|subject_id)`, not the omics ID and not the pedigree file's column name:

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics expression.csv `
  --id omics_id --individual-id subject_id `
  --formula "BMI ~ Sex + Age + <omics> + (1|subject_id)" `
  --pedigree pedigree.csv --pedigree-id member_id `
  --sire-id parent1_id --dam-id parent2_id `
  --out bmi-pedigree-expression.csv
```

Blank, `NA`, `0`, `.`, and `-9` parent values mean unknown. Pedigree
members need not have phenotypes; they can be needed to connect relatives.
Subject `9009` is absent from the example pedigree, so JLinAlg adds it as an
unrelated singleton and reports this on screen, in the log, and in the
manifest. Repeated phenotype rows with the same absent subject ID share one
singleton family and random effect.

Use `--pedigree-family-id` only if member IDs repeat between families. When
it is used, known phenotype IDs must be qualified as
`family_id:member_id`, for example `F10:1003`. Do not supply raw `1003`
for a known qualified member; that is a different ID and becomes a singleton.
The pedigree term uses Henderson's additive relationship precision. Another
term such as `(1|clinic)` remains an independent random intercept.

## 6. GLM and GLMM

Use `--family binomial` for a binary response, `poisson` for counts, and
`gamma` for positive continuous responses. Gaussian is the default.

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics expression.csv `
  --id omics_id --formula "case_status ~ Sex + Age + <omics>" `
  --family binomial --case-value case --control-value control `
  --out disease-expression.csv
```

Without a random effect this is a GLM. Add `(1|clinic)`, a GRM, or a matching
pedigree term for a first-order Laplace, glmer-like GLMM:

```powershell
java -jar jlinalg-<version>.jar --pheno phenotype.csv --omics expression.csv `
  --id omics_id --individual-id subject_id `
  --formula "case_status ~ Sex + Age + <omics> + (1|subject_id)" `
  --family binomial --case-value case --control-value control `
  --pedigree pedigree.csv --pedigree-id member_id `
  --sire-id parent1_id --dam-id parent2_id `
  --out disease-pedigree-expression.csv
```

`--link` overrides a supported family default when needed; omit it for the
canonical link. Genotype Laplace-GLMM scans are not currently exposed.

## 7. Cox survival analysis

The phenotype table needs a follow-up time and an event indicator. Use
`Surv(time,event)`; event is zero for censoring and one for an observed
event.

```powershell
java -jar jlinalg-<version>.jar --pheno survival.csv --id sample_id `
  --formula "Surv(followup,event) ~ Sex + Age" `
  --model cox --ties efron --out survival-cox.csv
```

Use `--ties breslow` only when that convention is required. Add a GRM and
`--individual-id` for Gaussian kinship frailty. Streamed `<omics>` Cox,
start-stop data, formula frailty, and pedigree-correlated Cox are Java-API
features rather than options in this general CLI.

## 8. Ridge, lasso, and elastic net

Penalized regression uses a numeric observation-by-column CSV/TSV:

```text
y,x1,x2,x3
2.1,0.2,1.1,0.0
3.4,0.9,0.7,1.0
```

```powershell
java -jar jlinalg-<version>.jar penalized-regression --input predictors.csv `
  --response y --predictors x1,x2,x3 --model ridge `
  --lambda-grid 10,3,1,0.3,0.1 --cv-folds 5 --out ridge.tsv

java -jar jlinalg-<version>.jar penalized-regression --input predictors.csv `
  --response y --predictors x1,x2,x3 --model lasso `
  --lambda-grid 1,0.3,0.1,0.03,0.01 --cv-folds 5 --out lasso.tsv

java -jar jlinalg-<version>.jar penalized-regression --input predictors.csv `
  --response y --predictors x1,x2,x3 --model elastic-net --alpha 0.5 `
  --lambda-grid 1,0.3,0.1,0.03,0.01 --cv-folds 5 --out elastic-net.tsv
```

Ridge fixes alpha to 0, lasso fixes it to 1, and elastic net uses the supplied
alpha between them. Predictors are standardized and an intercept is fitted by
default. Use `--no-standardize` or `--no-intercept` deliberately.
Coefficient rows do not carry classical p-values: selecting a penalty and
performing post-selection inference are different tasks.

## Next tutorials

- [CLI-only Mendelian randomization](cli-mr-tutorial.md)
- [CLI-only SuSiE fine mapping](cli-susie-tutorial.md)
- [CLI-only multi-signal colocalization](cli-colocalization-tutorial.md)
- [Full command reference](command-line.md)
