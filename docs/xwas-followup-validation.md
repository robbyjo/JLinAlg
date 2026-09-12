# xWAS follow-up workflows: validation and supported scope

Source-build additions, September 2026. These commands are not present in the
previously published v0.3.5 executable; build the current source.

Final local gate: `check javadoc executableJar` passed on 2026-09-12 with
772 discovered tests: **769 passed, three optional CHOLMOD skips, zero failures
or errors**. All 14 new workflow tests passed. Packaged-JAR smoke tests cover
all six commands, including shared LDSC-to-factor files and score training to
held-out application. The four web vignettes were inspected in the browser;
the static website link/structure checker also passed.

Numerical assertion tolerances are absolute: 1e-11 for LDSC S/V/intercepts,
1e-12 for molecular tests, 2e-6 for factor parameters and 1e-8 for factor
parameter covariance, and 1e-8 for ridge coefficients and held-out metrics.
The independent factor optimizer and Java score-certified solution agree
within those limits; no external package is treated as an infallible oracle.

| Workflow | CLI | Implemented estimator |
| --- | --- | --- |
| Genetic architecture | `ldsc` | One-step unpartitioned observed-scale LD-score regression, free intercepts, two weight updates, common final-weight block jackknife |
| Predicted molecular association | `twas`, `pwas` | Raw-dosage prediction-weight Z test with reference genotype SD/LD and optional full-rank joint model test |
| Joint genetic traits | `genomic-factor` | Full-WLS single-factor measurement model and GLS SNP tests conditional on fitted loadings |
| Prediction | `score-train`, `score-apply` | Gaussian ridge/elastic-net training/CV, portable original-scale weights, allele-aware application, held-out quantitative metrics or binary rank AUC |

## Independent numerical fixtures

`src/test/resources/xwas/reference.R` regenerates the checked-in fixtures and
`examples/xwas` inputs using base R. The script uses a fixed seed for synthetic
input creation; inferential methods themselves are deterministic. No real
cohort analysis is used as a validation oracle.

- LDSC: four traits, 1,200 variants, 30 shared blocks. Independent `lm.wfit`
  QR solves reproduce genetic S, intercepts and full jackknife sampling V.
  Regression weights follow the documented moment model. This validates the
  stated estimator, not every upstream LDSC mode or dataset calibration.
- TWAS/PWAS: direct matrix cross-products and normal tails, plus independent
  joint quadratic-form calculations. A deliberately swapped GWAS allele
  validates sign alignment. Tests exercise singular PSD LD, zero predicted
  variance, redundant models, incomplete coverage and invalid alleles.
- Genetic factor: base-R BFGS with numerically differentiated gradients fits
  a correlated sampling-covariance fixture; a separately evaluated numerical
  Jacobian supplies parameter covariance. Full-C SNP GLS estimates and
  heterogeneity are checked against direct R linear solves. The tests also
  pass actual LDSC S/V sidecars into the factor CLI.
- Scores: ridge estimates are checked against augmented-design QR in R,
  independently of coordinate descent. Model serialization/application and
  held-out RMSE, predictive R2 and calibration match direct R calculations.
  Tests verify deterministic CV, training-ID overlap rejection, allele swaps
  using 2-dosage, AUC ties and negative out-of-sample predictive R2.

Output-path collisions, label reordering, help dispatch and lifecycle logs are
covered by `XwasWorkflowTest`. All derived outputs are checked for existing
files before publication; command failures do not publish successful estimates.

Additional analytic tests exercise variable-N LDSC, negative heritability,
molecular weight-scale/orientation invariance, and genotype-to-molecular score
training followed directly by the TWAS command.

## Synthetic workload timings

Measured on the development Windows host with Java 25, a 1 GiB heap limit,
two warmups and five measured runs. Values are median in-process elapsed
times, excluding JVM startup, file parsing and writing. Inputs are generated
with seed 421; no real participant data were processed.

| Workload | Median |
| --- | --- |
| LDSC: 20,000 variants, four traits, 200 jackknife blocks | 19.250 ms |
| Molecular Z test: 100 reference variants, including LD validation | 0.752 ms |
| Genetic factor: four traits with supplied S/V | 0.545 ms |
| Ridge score: 2,000 samples, 50 features, one prespecified lambda | 4.685 ms |

Reproduce with `./gradlew.bat benchmarkXwasFollowup`. The checksum from the
measured workload is `14.68429368278705`. These small synthetic timings are
not a speed comparison against LDSC, MetaXcan, GenomicSEM or PRSice, and do not
predict end-to-end biobank or genome-wide score runtime.

## Reproduce

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/xwas/reference.R
./gradlew.bat test --tests org.jlinalg.cli.XwasWorkflowTest
./gradlew.bat check javadoc executableJar
```

In restricted Windows environments, use an accessible workspace temporary
directory for Java test processes and a writable Gradle user home. This changes
where JUnit creates temporary files, not numerical tolerances or test behavior.

## Interpretation and limits

The [LDSC vignette](vignettes/ldsc.md), [predicted-omics vignette](vignettes/predicted-omics.md),
[factor vignette](vignettes/genomic-factor.md) and [score vignette](vignettes/prediction-scores.md)
document schemas, equations, examples and exclusions. Important boundaries:

- Complete harmonized summaries and ancestry/scale-compatible reference data
  are required. Wide LDSC/factor inputs cannot detect an unreported per-trait
  allele reversal. Strand complements, genome lifts and ambiguity resolution
  are upstream work.
- LDSC does not include partitioned/two-step estimators, liability conversion,
  LD-score computation or summary munging. Correlations with nonpositive full
  or deleted heritabilities are explicitly unavailable. Estimated matrices
  are not made positive definite by automatic repair.
- Molecular tests require complete model coverage; joint tests reject
  rank-deficient model correlation. Native model-database formats and
  S-MultiXcan rank truncation are not included.
- Genomic factor fitting supports one factor with positive residual variances.
  SNP standard errors condition on loadings and do not propagate measurement
  uncertainty or SNP–measurement cross-covariance. Full GenomicSEM parity,
  robust DWLS, and multiple factors are not claimed.
- Gaussian score training does not implement LD-aware Bayesian PRS, logistic
  training or grouped CV. Binary evaluation is rank discrimination, not
  absolute-risk calibration. Independent cohort membership cannot be proven
  solely from sample ID checks, especially for imported external weights.
- Tables and small covariance matrices are read in memory. LDSC avoids a full
  variant rescan per jackknife block; prediction training remains a dense
  participant-by-feature workload, and molecular LD validation is cubic in
  locus variant count. This is not a genome-scale streaming PRS implementation.

## Method references

- [LDSC](https://github.com/bulik/ldsc): regression moments, weights and block jackknife.
- [MetaXcan](https://github.com/hakyimlab/MetaXcan): genetically predicted molecular association.
- [Genomic SEM](https://www.nature.com/articles/s41562-019-0566-x): genetic and sampling covariance approach.
- [PRSice](https://github.com/choishingwan/PRSice): polygenic workflow context.

Package names identify method references and intended interoperability, not
complete feature parity. The five additional analysis families requested for
later work are tracked in [TODO.md](../TODO.md).
