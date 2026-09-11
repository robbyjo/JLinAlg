# Command-line-only workflows

Every workflow here is invoked through the executable JAR; no Java API code or
dataframe runtime is required. Start a new general association
file/formula/pedigree with <code>--dry-run</code>. The commands and options on
this page are available in `jlinalg-0.3.5.jar`.

For a slower first-principles walkthrough, start with the
[progressive association tutorial](cli-association-tutorial.md). Separate
file-to-result tutorials cover [mediation](cli-mediation-tutorial.md),
[Mendelian randomization](cli-mr-tutorial.md),
[SuSiE](cli-susie-tutorial.md), and
[multi-signal colocalization](cli-colocalization-tutorial.md).

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
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "BMI ~ Sex + Age + <omics>" --out bmi-expression.csv --dry-run
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

## Phenotype-only OLS

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --id SampleName --formula "BMI ~ Sex + Age" --out bmi-coefficients.csv
~~~

Gaussian data with no random term, GRM, or pedigree selects OLS. Rows missing
the response or any model variable are omitted jointly.

## Numeric omics OLS

Numeric omics files are feature-by-sample: the first column contains feature
IDs and the remaining headers contain sample IDs.

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --omics-type expression --formula "BMI ~ Sex + Age + <omics>" --transform "<omics> = winsor_mad(k=4) | zscore()" --annot genes.tsv --annot-id gene_id --annot-cols symbol,chromosome --threads 8 --block-size auto --out bmi-expression.csv
~~~

Phenotype and omics IDs are intersected automatically. Missing phenotype rows
are omitted; missing values within an omics feature are mean-imputed after
transformation.

## Grouped numeric LMM

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "BMI ~ Sex + Age + <omics> + (1|site)" --variance-components auto --df satterth --out bmi-expression-site.tsv
~~~

For numeric omics, <code>auto</code> refits variance components by REML for
every feature. This is the exact lmer-like path, not P3D/EMMAX.

## Pedigree LMM with globally unique IDs

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --individual-id sabreid --formula "BMI ~ Sex + Age + <omics> + (1|sabreid)" --pedigree pedigree.csv --pedigree-id sabreid --sire-id fid --dam-id mid --out bmi-expression-pedigree.csv
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
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --id SampleName --individual-id pedigree_key --formula "BMI ~ Age + (1|pedigree_key)" --pedigree pedigree.csv --pedigree-family-id family --pedigree-id member --sire-id sire --dam-id dam --out bmi-family-pedigree.tsv
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
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --individual-id sabreid --formula "BMI ~ Sex + Age + <omics> + (1|sabreid) + (1|batch)" --pedigree pedigree.csv --pedigree-id sabreid --sire-id fid --dam-id mid --out bmi-pedigree-batch.tsv
~~~

Only the term matching <code>--individual-id</code> receives pedigree
precision. <code>(1|batch)</code> is an independent random intercept.

## GRM-adjusted numeric scan

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --individual-id IID --formula "BMI ~ Sex + Age + <omics>" --grm cohort --out bmi-expression-grm.tsv
~~~

<code>--grm</code> accepts a labeled square CSV/TSV or a GCTA prefix backed by
<code>cohort.grm.bin</code> and <code>cohort.grm.id</code>. Numeric features get
exact per-feature REML refits. A GRM and pedigree cannot be combined.

## Binary, count, and positive outcomes

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "case_status ~ Sex + Age + <omics>" --family binomial --case-value case --control-value control --out disease-expression.tsv
~~~

Add <code>(1|site)</code>, a GRM, or a matching pedigree term for a
first-order Laplace GLMM with per-feature glmer-like refits. Use
<code>--family poisson</code> for counts and <code>--family gamma</code> for
positive continuous outcomes. Genotype Laplace GLMM scans are not exposed.

## Genotype GWAS and P3D LMM

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --omics cohort.vcf.gz --id IID --formula "BMI ~ Sex + Age + PC1 + PC2 + <omics>" --min-maf 0.01 --min-mac 20 --max-marker-missing 0.02 --out bmi-gwas.tsv
~~~

BCF and biallelic layout-2 BGEN are accepted; add
<code>--sample-file cohort.sample</code> when BGEN IDs are external. Adding a
GRM selects genotype LMM null-model P3D/EMMAX. Pedigree genotype scans and
per-marker variance-component refits are not yet available.

## Cox regression and GRM frailty

~~~powershell
java -jar jlinalg-0.3.5.jar --pheno phenotype.csv --id SampleName --formula "Surv(followup,event) ~ Sex + Age" --ties efron --out survival.tsv
~~~

Add <code>--grm cohort --individual-id IID</code> for Gaussian kinship frailty.
Formula frailty terms and pedigree Cox are available through the Java API but
not yet through this general CLI.

## Specialized CLI subcommands

~~~powershell
java -jar jlinalg-0.3.5.jar beta-regression --input proportions.tsv --response proportion --mean dose,age --precision batch --out beta.tsv
java -jar jlinalg-0.3.5.jar penalized-regression --input continuous.tsv --response y --predictors x1,x2,x3 --model elastic-net --alpha 0.5 --lambda-grid 1,0.3,0.1,0.03 --cv-folds 5 --out elastic-net.tsv
java -jar jlinalg-<version>.jar mediation --input mediation.tsv --outcome y --treatment x --mediator m --covariates age --out mediation-effects.tsv
java -jar jlinalg-<version>.jar susie --summary locus.tsv --ld locus-ld.tsv --sample-size 10000 --out locus-susie.tsv
java -jar jlinalg-<version>.jar coloc --trait1 trait1-susie.tsv.effects.tsv --trait2 trait2-susie.tsv.effects.tsv --out coloc.tsv
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
