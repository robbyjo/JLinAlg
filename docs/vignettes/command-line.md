# Command-line-only workflows

JLinAlg v0.3.6 additionally provides [LDSC](ldsc.md),
[genetically predicted TWAS/PWAS](predicted-omics.md),
[genetic factors and common-factor GWAS](genomic-factor.md), and
[prediction score training/application](prediction-scores.md).
Use `ldsc`, `twas`, `pwas`, `genomic-factor`, `score-train`, or `score-apply`
with `--help` for their independent table-based schemas.

Every workflow here is invoked through the executable JAR; no Java API code or
dataframe runtime is required. Start a new general association
file/formula/pedigree with <code>--dry-run</code>. The complete documented CLI is
available in `jlinalg-0.3.6.jar` and later.

For a slower first-principles walkthrough, start with the
[progressive association tutorial](cli-association-tutorial.md). Separate
file-to-result tutorials cover [mediation](cli-mediation-tutorial.md),
[Mendelian randomization](cli-mr-tutorial.md),
[SuSiE](cli-susie-tutorial.md), and
[multi-signal colocalization](cli-colocalization-tutorial.md).

JLinAlg v0.3.6 also provides `meta-analysis` and `meta-regression`
for separate cohort summary files. See the [meta-analysis tutorial](cli-meta-analysis-tutorial.md)
for cohort directions, missing values, and minimum-cohort controls.

The release also provides [`grm` construction](../grm-cli.md) from genotype
files and [rare-variant meta-analysis](../rare-variant-meta-analysis.md).

Latent-factor and batch workflows use dedicated global-matrix commands rather
than the row-wise `--transform` pipeline. See the
[latent-confounder vignette](latent-confounders-and-batch.md) for PCA, SVA,
AutoSVA, PEER, and ComBat commands. The recently added probit, prediction,
IV/2SLS, conditional-score, and ARIMA-smoothing surfaces are available as
`--family probit`, `glm-predict`, `iv-regression`, `conditional-score`, and
`arima-regression`.

## Understand the three ID layers

| Option | Meaning |
| --- | --- |
| <code>--id SampleName</code> | Phenotype observation/sample ID. Omics sample IDs are intersected and reordered against it. |
| <code>--individual-id sabreid</code> | Phenotype column mapping observations to GRM or pedigree members. It may repeat and defaults to <code>--id</code>. |
| <code>--pedigree-id sabreid</code> | Individual column inside the pedigree file; sire and dam columns refer to these values. |

The pedigree term must use the phenotype matching-column name:
<code>--individual-id sabreid</code> requires <code>(1|sabreid)</code>.
A different term such as <code>(1|batch)</code> stays an ordinary independent
random effect.

After sample alignment and complete-case omission, phenotype individual IDs
absent from the pedigree file are added as unrelated, noninbred singleton
families. Repeated observations sharing one missing individual ID map to the
same singleton founder; different missing IDs remain unrelated. The screen,
log, and manifest report singleton-family and singleton-observation counts.
This retains observations without inventing parentage. At least one aligned
phenotype ID must resolve to a source pedigree member; a zero-match run fails
rather than replacing the requested pedigree with an all-singleton model.

## Inspect routing without fitting

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "BMI ~ Sex + Age + <omics>" --out bmi-expression.csv --dry-run
~~~

For the general command, <code>--dry-run</code> reads and aligns inputs,
selects complete phenotype rows, compiles the formula, loads a requested GRM
or pedigree, validates the transform, and prints model routing, output format,
backend policy, and (for omics) block/worker sizing. It then stops before
model fitting or feature scanning, leaving only the normal log unless
<code>--no-log</code> is set. <code>--explain</code> prints the same plan and
continues with the fit. The option is singular: <code>--explain</code>.
Specialized subcommands have their own <code>--help</code> and do not accept
these two general-command switches.
Logs default to <code>OUT.log</code>; manifests default to
<code>OUT.manifest.json</code>. A <code>.csv</code> output is comma-separated;
<code>.tsv</code> and other suffixes are tab-separated.

## Run timestamps and elapsed time

Source builds now include these fields in CLI logs:

```text
started=2026-09-12T14:00:00Z
finished=2026-09-12T16:03:04.125Z
elapsed_ms=7384125
elapsed=2 h 3 min 4.125 s
status=complete
```

Timestamps are UTC; elapsed time uses a monotonic clock. Durations include
milliseconds, and days for runs lasting at least 24 hours. The start record is
flushed before computation. Once execution has begun, handled failures retain
the log with a finish timestamp, elapsed time, and `status=failed`. Each appended
run has its own run ID and timing block.

The general association CLI retains its timestamped lab-notebook format;
meta-analysis, meta-regression, mediation, SuSiE, and colocalization retain
their command-specific metadata plus these timing fields. Existing explicit
log paths and the general CLI's `--no-log` option continue to work.

`ld-db`, `mr-instruments`, `clump`, `mr-xwas`, `mr-estimate`, `mr-mvmr`,
`mr-multivariate`,
`beta-regression`, and `penalized-regression` now also write run logs when
invoked through the executable. Their defaults are `OUT.log` when an
`--out`/`--output` is supplied, or a unique `logs/jlinalg-COMMAND-*.log`
otherwise. Use `--log FILE` to choose a path or `--no-log` to opt out.
The log location is printed to stderr, preserving machine-readable stdout.
Their existing logs are protected unless `--overwrite` is requested, in which
case another timing block is appended. Help/version requests do not create logs.

## Phenotype-only OLS

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --id SampleName --formula "BMI ~ Sex + Age" --out bmi-coefficients.csv
~~~

Gaussian data with no random term, GRM, or pedigree selects OLS. Rows missing
the response or any model variable are omitted jointly.

## Numeric omics OLS

Numeric omics files are feature-by-sample: the first column contains feature
IDs and the remaining headers contain sample IDs.

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --omics-type expression --formula "BMI ~ Sex + Age + <omics>" --transform "<omics> = winsor_mad(k=4) | zscore()" --annot genes.tsv --annot-id gene_id --annot-cols symbol,chromosome --threads 8 --block-size auto --out bmi-expression.csv
~~~

Phenotype and omics IDs are intersected automatically. Missing phenotype rows
are omitted; missing values within an omics feature are mean-imputed after
transformation.

## Grouped numeric LMM

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "BMI ~ Sex + Age + <omics> + (1|site)" --variance-components auto --df satterth --out bmi-expression-site.tsv
~~~

For numeric omics, <code>auto</code> refits variance components by REML for
every feature. This is the exact lmer-like path, not P3D/EMMAX.

## Pedigree LMM with globally unique IDs

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --individual-id sabreid --formula "BMI ~ Sex + Age + <omics> + (1|sabreid)" --pedigree pedigree.csv --pedigree-id sabreid --sire-id fid --dam-id mid --out bmi-expression-pedigree.csv
~~~

The <code>(1|sabreid)</code> term uses Henderson's sparse additive relationship
inverse. Pedigree members without phenotypes remain in the ancestry graph.
Phenotype IDs missing from the pedigree become reported singleton families;
repeated rows with the same ID share one founder/random effect.
Blank, <code>NA</code>, <code>0</code>, <code>.</code>, and <code>-9</code>
parent IDs mean unknown.

## Family-disambiguated pedigree IDs

Use <code>--pedigree-family-id family</code> to disambiguate individual IDs that
repeat between families. A globally unique member still matches its raw ID,
and the exact <code>family:individual</code> form is accepted as an alias. If a
raw member ID occurs in more than one family, however, an unqualified
phenotype value is ambiguous and the run fails with a request to use the exact
qualified key:

~~~text
SampleName,pedigree_key,BMI,Age
S1,F01:1001,24.1,51
S2,F01:1002,26.7,49
S3,SINGLETON_3,22.9,61
~~~

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --id SampleName --individual-id pedigree_key --formula "BMI ~ Age + (1|pedigree_key)" --pedigree pedigree.csv --pedigree-family-id family --pedigree-id member --sire-id sire --dam-id dam --out bmi-family-pedigree.tsv
~~~

<code>F01:1001</code> maps to the qualified pedigree member.
<code>SINGLETON_3</code>, which has no family/pedigree row, becomes a singleton.
Raw <code>1001</code> also maps when it occurs only once in the pedigree. If both
<code>F01:1001</code> and <code>F02:1001</code> exist, use the qualified value.

Family labels are namespaces, not relationship groups. Parent-child links
construct the ancestry graph. A unique raw parent reference can therefore
connect a child to a parent whose row has another family label; conversely,
sharing a family label alone does not make two founders related.

## Pedigree plus an ordinary grouping effect

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --individual-id sabreid --formula "BMI ~ Sex + Age + <omics> + (1|sabreid) + (1|batch)" --pedigree pedigree.csv --pedigree-id sabreid --sire-id fid --dam-id mid --out bmi-pedigree-batch.tsv
~~~

Only the term matching <code>--individual-id</code> receives pedigree
precision. <code>(1|batch)</code> is an independent random intercept.

## GRM-adjusted numeric scan

To construct the input matrix first, use the source-build `grm` command:

```shell
java -jar build/cli/jlinalg-0.3.6.jar grm --genotypes cohort.vcf.gz --maf 0.01 --call-rate 0.95 --out cohort-grm.tsv
```

Pass `--grm cohort-grm.tsv` to the analysis command. See the
[construction tutorial](../grm-cli.md) for a runnable small example, input
formats, normalization, and dense-memory requirements.

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --individual-id IID --formula "BMI ~ Sex + Age + <omics>" --grm cohort --out bmi-expression-grm.tsv
~~~

<code>--grm</code> accepts a labeled square CSV/TSV or a GCTA prefix backed by
<code>cohort.grm.bin</code> and <code>cohort.grm.id</code>. Numeric features get
exact per-feature REML refits. A GRM and pedigree cannot be combined.

## Binary, count, and positive outcomes

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "case_status ~ Sex + Age + <omics>" --family binomial --case-value case --control-value control --out disease-expression.tsv

java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "case_status ~ Sex + Age + <omics>" --family probit --case-value case --control-value control --out disease-expression-probit.tsv
~~~

Add <code>(1|site)</code>, a GRM, or a matching pedigree term for a
first-order Laplace GLMM with per-feature glmer-like refits. The source build's
<code>--family probit</code> uses a stable binomial normal-CDF link; use
<code>--family poisson</code> for counts and <code>--family gamma</code> for
positive continuous outcomes. Genotype Laplace GLMM scans are not exposed.

## Genotype GWAS and P3D LMM

The current source build adds `--conditional-gwas-summary` with
`--score-genome-build` for Gaussian, logistic, probit, Poisson and model-based Cox
genotype scans. It writes aggregate scores, bounded covariance blocks and
coding/model metadata. Add `--condition-on` to refit the cohort null with lead
variants. See the [conditional-GWAS workflow](../conditional-gwas-summary.md)
for complete commands, column definitions and calibration/coverage limits.
Complete blocks can then be imported by the
[summary-only conditioning workflow](conditional-score-conditioning.md).

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --omics cohort.vcf.gz --id IID --formula "BMI ~ Sex + Age + PC1 + PC2 + <omics>" --min-maf 0.01 --min-mac 20 --max-marker-missing 0.02 --out bmi-gwas.tsv
~~~

BCF and biallelic layout-2 BGEN are accepted; add
<code>--sample-file cohort.sample</code> when BGEN IDs are external. Adding a
GRM selects genotype LMM null-model P3D/EMMAX. Pedigree genotype scans and
per-marker variance-component refits are not yet available.

## Cox regression and GRM frailty

~~~powershell
java -jar jlinalg-0.3.6.jar --pheno phenotype.csv --id SampleName --formula "Surv(followup,event) ~ Sex + Age" --ties efron --out survival.tsv
~~~

Add <code>--grm cohort --individual-id IID</code> for Gaussian kinship frailty.
Formula frailty terms and pedigree Cox are available through the Java API but
not yet through this general CLI.

## Specialized CLI subcommands

~~~powershell
java -jar jlinalg-0.3.6.jar beta-regression --input proportions.tsv --response proportion --mean dose,age --precision batch --out beta.tsv
java -jar jlinalg-0.3.6.jar penalized-regression --input continuous.tsv --response y --predictors x1,x2,x3 --model elastic-net --alpha 0.5 --lambda-grid 1,0.3,0.1,0.03 --cv-folds 5 --out elastic-net.tsv
java -jar jlinalg-<version>.jar mediation --input mediation.tsv --outcome y --treatment x --mediator m --covariates age --out mediation-effects.tsv
java -jar jlinalg-<version>.jar susie --summary locus.tsv --ld locus-ld.tsv --sample-size 10000 --out locus-susie.tsv
java -jar jlinalg-<version>.jar coloc --trait1 trait1-susie.tsv.effects.tsv --trait2 trait2-susie.tsv.effects.tsv --out coloc.tsv
java -jar jlinalg-<version>.jar confounders --method sva --omics expression.tsv --pheno phenotype.tsv --id IID --full-design case,age,sex --null-design age,sex --factors auto --out expression-sva
java -jar jlinalg-<version>.jar batch-adjust --method combat --omics expression.tsv --pheno phenotype.tsv --id IID --batch plate --preserve case,age,sex --out expression-combat
java -jar jlinalg-<version>.jar glm-predict --input model.tsv --response y --predictors exposure,age --family probit --estimand expected --out predictions.tsv
java -jar jlinalg-<version>.jar iv-regression --input iv.tsv --response y --exogenous age --endogenous exposure --instruments z1,z2 --out iv-fit
java -jar jlinalg-<version>.jar arima-regression --input series.tsv --response y --predictors time --order 1,0,1 --smooth --out arima-fit
~~~

MR preparation, LD installation/clumping, MR estimation, and parallel xWAS MR
are CLI-native. See the dedicated
[mediation](cli-mediation-tutorial.md), [MR](cli-mr-tutorial.md),
[SuSiE](cli-susie-tutorial.md), and
[colocalization](cli-colocalization-tutorial.md) tutorials and the parallel
[xWAS MR guide](xwas-mr-pipeline.md) for their schemas.

## Operational checklist

- Start a general association command with <code>--dry-run</code>; use the
  singular <code>--explain</code> to print routing and still fit. These
  switches do not apply to specialized subcommands.
- Automatic block sizing uses both JVM heap headroom and the requested thread
  count. When memory permits it queues at least two complete worker waves so
  work stealing can absorb slow chunks at block boundaries. Low heap headroom
  uses the largest complete wave that fits or safely limits worker capacity;
  explicit <code>--block-size N</code> remains available.
- Inspect aligned, omitted, and pedigree-singleton counts before interpreting
  results.
- Keep <code>--id</code> for observation/omics alignment and
  <code>--individual-id</code> for pedigree or GRM membership.
- Use <code>--overwrite</code> deliberately. Verified partial-output resume is
  not available; use a new output path after interruption.
- Preserve the generated log and manifest with the result table.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [J. A. Nelder and R. W. M. Wedderburn (1972) — Generalized linear models](../CITATIONS.md#nelder-wedderburn-1972)
- [H. D. Patterson and Robin Thompson (1971) — Recovery of inter-block information when block sizes are unequal](../CITATIONS.md#patterson-thompson-1971)
- [D. R. Cox (1972) — Regression models and life-tables](../CITATIONS.md#cox-1972)
- [Norman E. Breslow and David G. Clayton (1993) — Approximate inference in generalized linear mixed models](../CITATIONS.md#breslow-clayton-1993)
- [Hyun Min Kang et al. (2010) — Variance component model to account for sample structure in genome-wide association studies](../CITATIONS.md#kang-emmax-2010) — [PMID: 20208533](https://pubmed.ncbi.nlm.nih.gov/20208533/) · [PMCID: PMC3092069](https://pmc.ncbi.nlm.nih.gov/articles/PMC3092069/)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
