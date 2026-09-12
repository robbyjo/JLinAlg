# Synthetic xWAS examples

These tables contain synthetic data only. No participant data or licensed
prediction models are included. `src/test/resources/xwas/reference.R`
recreates the inputs and independent base-R numerical references.

Build the current source with `./gradlew.bat executableJar`, then follow:

- [LDSC](../../docs/vignettes/ldsc.md): `ldsc.tsv` contains four harmonized
  synthetic traits, reference LD scores, weighting LD scores and sample sizes.
- [Predicted TWAS/PWAS](../../docs/vignettes/predicted-omics.md): `gwas.tsv`,
  `weights.tsv`, `reference.tsv` and `ld.tsv` define two molecular predictors.
- [Genetic factor](../../docs/vignettes/genomic-factor.md): `factor-S.tsv`,
  `factor-V.tsv`, `factor-gwas.tsv` and `sampling-correlation.tsv` separate
  genetic covariance, moment uncertainty, and SNP sampling overlap.
- [Scores](../../docs/vignettes/prediction-scores.md): `train.tsv` and `test.tsv`
  are disjoint quantitative cohorts. The three `polygenic-weights.tsv`,
  `target-alleles.tsv`, and `target-dosages.tsv` tables demonstrate an imported
  score, an allele reversal and binary rank-AUC evaluation.

Outputs should go into fresh paths under `build/xwas/`. Tutorial p-values and
performance do not establish calibration or predictive utility in real cohorts.
