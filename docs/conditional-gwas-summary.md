# Cohort-side conditional GWAS summaries

The general genotype CLI can export efficient scores and their covariance for
unrelated-sample Gaussian, binary logistic or probit, Poisson-log, and model-based Cox
analyses. Cohorts retain the phenotype and genotype records. The shared files
contain aggregate statistics and model metadata. A separate
`conditional-score` command now validates and conditions these summaries; see
the [summary-only conditioning vignette](vignettes/conditional-score-conditioning.md).

This is available in the current source build; build with
`./gradlew.bat executableJar`. It does not change the existing `mr-estimate`
Gaussian summary model or make ordinary marginal Wald statistics sufficient to
reconstruct a logistic or Cox likelihood.

## Export from a cohort

```powershell
java -jar build/cli/jlinalg-0.3.6.jar `
  --omics cohort.vcf.gz --pheno phenotype.tsv --id IID `
  --formula "case_status ~ age + sex + PC1 + PC2 + <omics>" `
  --family binomial --case-value case --control-value control `
  --conditional-gwas-summary --score-genome-build GRCh38 `
  --score-block-size 64 --out binary-gwas.tsv
```

Use `--family poisson` for counts, with `offset(log_exposure)` in the formula
when appropriate. Gaussian OLS and Gaussian GLM exports are also supported.
Cox genotype scans now accept the ordinary survival formula:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar `
  --omics cohort.vcf.gz --pheno phenotype.tsv --id IID `
  --formula "Surv(entry,time,event) ~ age + sex + PC1 + <omics>" `
  --ties efron --conditional-gwas-summary --score-genome-build GRCh38 `
  --out cox-gwas.tsv
```

`Surv(time,event)` supports right censoring; the three-column form also supports
left truncation. Event values must be 0/1. Efron and Breslow ties are supported.
The current CLI requires one observation per sample and at least one varying
null covariate or conditioning variant. The compiler uses ordinary reference
contrasts for categorical covariates, then removes the Cox intercept.
Strata, repeated subject intervals, formula weights, frailty, GRM, pedigree,
and cluster/relatedness score covariance are outside this Cox export path.

The switch adds three linked outputs to the normal run:

| File | Content |
| --- | --- |
| `binary-gwas.tsv` | Existing columns plus the score fields below |
| `binary-gwas.tsv.score-cov.tsv` | Upper triangle, including the diagonal, of each complete score-covariance block |
| `binary-gwas.tsv.score-manifest.json` | Schema, model, formula, covariate and conditioning identities, genome build, coding, dispersion, convergence, calibration and coverage |

The ordinary run log and manifest are retained. Non-Gaussian genotype scans
print and log a notice that standard output is insufficient for MR conditional
GWAS. With the switch, the notice explains the local-score and calibration
limits. Gaussian scans and phenotype-only/numeric-omics regressions do not
receive that warning. `--no-log` suppresses the log file, retaining the screen
notice and the independently sufficient coding metadata in the score manifest.

## Extra columns and their units

| Columns | Meaning |
| --- | --- |
| `score_variant_key` | `CHR:POS:REF:ALT`, interpreted with the manifest's genome build; ALT is always the effect allele |
| `score_u`, `score_variance` | Unstandardized, nuisance-adjusted efficient score and model-based score variance at the fitted null |
| `n_analyzed` | Fixed number of aligned, complete phenotype/covariate rows, including mean-imputed genotypes |
| `n_cases`, `n_controls` | Binary analysis-sample counts, including imputed genotypes |
| `effect_ac_cases`, `effect_ac_controls` | ALT dosage sums among called genotypes, before imputation; fractional dosages are allowed |
| `n_called_cases`, `n_called_controls` | Called sample denominators for those allele sums |
| `n_events` | Number of observed events for Cox |
| `beta_score`, `se_score` | One-step `U/V` and `1/sqrt(V)`; these are not refitted alternative-model maximum-likelihood estimates |
| `p_score_normal` | Two-sided standard-normal tail of `U/sqrt(V)` |
| `p_score_calibrated` | Empty: SPA or another rare-event calibration has not been implemented |
| `calibration_method`, `calibration_status` | `normal-score`, `normal_only` |
| `null_model_id`, `conditioning_set_id`, `score_covariance_block` | Join keys identifying the exact fitted null, its conditioning set and matrix coverage |

Inapplicable values are empty. Filtered or nonestimable rows retain their
ordinary status/reason and have empty score fields; they have no covariance
entries. Binary export requires individual 0/1 responses after case/control
mapping and unit trial weights; grouped proportions require a different schema.
The ordinary `p_value`, estimates and `fdr_bh` are unchanged by enabling the
switch. In particular, `fdr_bh` still adjusts ordinary `p_value`, and the GLM
scanner's existing approximate t tail differs from the added normal score tail.

For a canonical GLM, the export residualizes genotypes against the null
covariates in the working-weight metric. Both U and V include the fitted null
dispersion scale; binomial and Poisson dispersion is one. Gaussian dispersion
is the fitted null estimate, held fixed for score calculations. Cox uses the
partial-likelihood score and information, projecting out nuisance covariates
with the Schur complement, including any small residual nuisance score.
Offsets and the actual complete-case sample enter these calculations.

## Refit a specified conditioning set locally

Add the lead variants to the cohort run:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar `
  --omics cohort.vcf.gz --pheno phenotype.tsv --id IID `
  --formula "case_status ~ age + sex + PC1 + PC2 + <omics>" `
  --family binomial --conditional-gwas-summary --score-genome-build GRCh38 `
  --condition-on rs123,1:456789:A:G --out binary-conditioned.tsv
```

This performs a source prepass, resolves each requested ID or allele-specific
key, and mean-imputes its dosage on the actual analysis sample. It then appends
the lead dosages to the covariate design and refits the null likelihood before
testing the remaining variants. Conditioning variants are not retested. Missing
IDs, ambiguous IDs, duplicate aliases, rank-deficient null designs and
nonconverged nulls fail explicitly. Up to 512 conditioning variants are accepted.

The resulting target score tests are conditional on that fitted null. They
remain asymptotic score tests. Changing the conditioning set requires another
cohort-side refit for likelihood-based nonlinear conditioning. One exported
score vector and Hessian describe local curvature, not the full likelihood
away from the fitted null. The `conditional-score` importer below can instead
apply that local Schur complement to an already exported complete score block.
The existing `mr-estimate --condition-on` remains a separate Gaussian beta/SE
plus LD approximation with an unchanged interface; it does not import this
schema. These conditional SNP association p-values are not causal MR p-values.

## Import and condition a complete score block

List one or more independent cohort export trios in a tab-separated manifest:

```text
cohort	summary	covariance	manifest
cohort_a	cohort-a.tsv	cohort-a.tsv.score-cov.tsv	cohort-a.tsv.score-manifest.json
cohort_b	cohort-b.tsv	cohort-b.tsv.score-cov.tsv	cohort-b.tsv.score-manifest.json
```

Then request output-oriented allele keys explicitly:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar conditional-score `
  --cohorts cohorts.tsv `
  --targets 1:456789:A:G,1:456950:C:T `
  --condition-on 1:455100:G:A --out locus-conditional.tsv
```

The importer permits exact forward-strand matches. It aligns REF/ALT swaps only
for Cox or explicit-intercept null models, where dosage translation is removed
by the null score; otherwise it requires the exported orientation. The same
rule applies when comparing fitted conditioning sets across cohorts. It checks
the completed schema/manifest joins, compatible null contracts, full
selected-block upper triangles, covariance diagonals, and conditioning rank.
It applies the Schur complement within each cohort before summing independent
cohort scores. Output metadata labels this `local_schur_one_step` inference and
states that no nonlinear cohort refit or cohort-overlap correction was done.
Targets and conditions in different blocks fail because the missing covariance
is unknown. See the vignette for output columns, API use, and failure handling.

## Covariance coverage and operational contract

Variants must be chromosome/position sorted, with unique allele-specific keys.
Use a consistent chromosome naming convention, normalized biallelic REF/ALT,
and additive ALT dosage in [0,2]. The exporter checks key syntax/order and
dosage bounds; it does not perform genome liftover, reference-FASTA verification
or variant normalization. It rejects out-of-range/infinite dosages and imputes NaN
genotypes with the analysis-sample mean, matching the scan.

Each covariance block contains at most `--score-block-size` successfully tested
variants, defaults to 64, and cannot exceed 512 or cross chromosomes. This is
independent of the scan's `--block-size`. Storage for retained genotypes is
O(n b), plus O(b²) covariance and duplicate-detection keys at the current
genomic position, beyond the ordinary model/scan workspace.
Genomic distance does not determine the block boundary. Covariance between
different blocks is **unavailable, not zero**; do not substitute ordinary LD
for it. A downstream conditional calculation needs every requested target,
condition and cross-covariance entry under the same null. Choose a block size
that covers the desired region, restrict the input region, or use the local
refit route when a conditioning lead is outside the block.

All score controls require the switch; it also requires `--score-genome-build`.
Mixed, quasi-likelihood and other outcome models need dedicated exporters and
are rejected by this switch. Existing Gaussian REML `rare-score --grm` remains
a separate workflow. Numerically clamped GLM score weights are rejected.

Outputs are preflighted against input and output-path aliases. Existing files
require `--overwrite`; `--resume` is not supported for this linked export.
The score manifest is published last and is the export completion marker.
An interrupted export may leave partial files; do not consume them or a
summary/covariance pair without its completed, matching manifest. On overwrite,
an old score success manifest is removed before writing replacement results.
Manifest variant/covariate identities are JSON arrays; other values follow the
existing manifest convention of strings or null.

## Independent validation and reproduction

`ConditionalGwasExportTest` compares all score and covariance entries with
checked-in R 4.6.1 `glm.fit` and `survival::coxph` results. Four variants and 179
complete rows cover missing genotypes, an omitted covariate row, offsets,
Gaussian dispersion, binary allele/sample counts, Poisson, both Cox tie rules,
right censoring, left truncation, and independently refitted conditioning nulls.
The score/covariance absolute tolerance is 2e-6. Operational tests cover default
output parity, notices, no-log coding, block boundaries, filtered rows,
conditioning aliases, unsupported models, and failed-overwrite completion.
`ConditionalScoreInferenceTest` checks allele-swap alignment, cohort-wise Schur
conditioning, pooling, compatibility and rank failures. `ConditionalScoreCliTest`
round-trips the exported fixture, matches a hand Schur complement, and verifies
that incomplete or cross-block covariance is rejected rather than zero-filled.

```powershell
# Optional fixture regeneration; R with survival must already be installed.
Rscript src/test/resources/r-reference/generate-conditional-score-export.R
./gradlew.bat test --tests org.jlinalg.cli.ConditionalGwasExportTest
./gradlew.bat test --tests org.jlinalg.settest.ConditionalScoreInferenceTest
./gradlew.bat test --tests org.jlinalg.cli.ConditionalScoreCliTest
./gradlew.bat benchmarkConditionalScoreExport
./gradlew.bat check benchmarkClasses javadoc executableJar
```

The benchmark repeats the four fixture genotypes at 64 distinct variant keys
and verifies all 64 scores and 2,080 upper-triangle entries against R. It times
complete CPU CLI runs, including input/output and FDR, with three warmups and
seven measured runs. This synthetic workload measures export overhead; it is
not a genome-scale throughput or rare-case calibration study.

Local Windows CPU results, 2026-09-12, with the benchmark run separately from
the test suite:

| Model, 179 rows / 64 variants | Standard median ms | Export median ms | Maximum absolute covariance error vs R |
| --- | ---: | ---: | ---: |
| Logistic | 18.4332 | 26.1667 | 1.00e-12 |
| Poisson | 15.5655 | 22.9674 | 1.99e-13 |
| Gaussian | 15.7001 | 21.7423 | 2.13e-14 |
| Cox, Efron | 16.4182 | 25.4181 | 1.92e-13 |
| Cox, Breslow | 17.1476 | 23.2759 | 2.90e-8 |

The full gate ran 727 tests: 724 passed, three optional native CHOLMOD tests
were skipped, and none failed. Website checks, benchmark compilation, Javadoc
and executable-JAR generation also passed.
