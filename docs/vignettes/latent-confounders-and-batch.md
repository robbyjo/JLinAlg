# Latent confounders and known batch effects

JLinAlg provides four sample-factor estimators and a separate known-batch
adjuster. The matrix orientation is always **features by samples**: the first
column contains feature IDs and the remaining header fields are sample IDs.
Measurements must already be normalized and approximately continuous.

## Choose the method deliberately

| Method | Use it for | Required information |
|---|---|---|
| PCA | A fast unsupervised baseline | A fixed factor count |
| SVA | Unknown unwanted variation while protecting a biological design | Full and null designs; factor count or BE selection |
| AutoSVA | Roby Joehanes' iteratively reweighted SVA and preserved-signal factor search | A protected design; fixed or automatic factor count |
| PEER | Probabilistic factor analysis with optional observed covariates | A fixed maximum factor count |
| ComBat | Location/scale adjustment for a known categorical batch | Batch labels and optional covariates to preserve |

PCA, SVA, AutoSVA, and PEER return sample covariates. ComBat returns an adjusted
feature matrix. They are not interchangeable meanings of "batch correction."

## PCA

```powershell
java -jar jlinalg-0.3.5.jar confounders `
  --method pca `
  --omics expression.tsv `
  --factors 10 `
  --out results/expression-pca `
  --write-adjusted
```

The outputs are:

- `expression-pca.factors.tsv`: sample by PC;
- `expression-pca.loadings.tsv`: feature by PC;
- `expression-pca.adjusted.tsv`: optional residual matrix;
- `expression-pca.manifest.tsv`: dimensions, settings, backend and variance
  fractions.

Use `--scale` only when equal feature variances are scientifically appropriate.
With `--pheno`, `--id`, and `--full-design age,sex`, PCA is fitted after those
known covariates are projected out. Sample IDs are intersected in omics order;
the aligned count is printed.

Tall matrices use a sample-Gram decomposition through the selected JDistlib
backend. This avoids allocating a feature-by-feature covariance matrix. Factor
signs are made deterministic by requiring the largest absolute loading to be
positive.

## Standard SVA

The full design contains the variables whose effects must be retained. The
null design contains adjustment variables but omits the biological term being
tested:

```powershell
java -jar jlinalg-0.3.5.jar confounders `
  --method sva `
  --omics methylation.tsv `
  --pheno phenotype.tsv `
  --id IID `
  --full-design case,age,sex `
  --null-design age,sex `
  --factors auto `
  --permutations 20 `
  --iterations 5 `
  --seed 7201 `
  --out results/methylation-sva
```

`--factors auto` follows `sva::num.sv(method="be")`: it residualizes by the
full model, independently permutes each feature across samples, and applies the
ordered 0.10 permutation threshold to residual singular-value fractions.
Numeric `--factors K` skips selection.

The iterative estimator follows `sva::irwsva.build`: nested-model F tests,
empirical local-FDR control probabilities, feature reweighting, and repeated
sample-factor decomposition. `feature-weights.tsv` makes the final empirical
control weights auditable.

## AutoSVA

```powershell
java -jar jlinalg-0.3.5.jar confounders `
  --method autosva `
  --omics methylation.tsv `
  --pheno phenotype.tsv `
  --id IID `
  --full-design case,age,sex `
  --null-design age,sex `
  --factors auto `
  --max-factors 30 `
  --increment 5 `
  --iterations 20 `
  --mrse 0.001 `
  --out results/methylation-autosva
```

AutoSVA is a named adaptation, not an alias for Bioconductor SVA. Its feature
weight is the probability of latent heterogeneity times the local-FDR
probability of not being driven by the protected model. Iteration stops at the
configured weight MRSE. Automatic selection evaluates the preserved-design F
ratio after residualizing each candidate factor count. The complete search,
including failures, is written to `factor-selection.tsv`.

The Java port deliberately makes degenerate-statistic handling deterministic.
It clamps invalid zero F statistics to a small positive value rather than
injecting random noise. It also interprets the source stopping expression as
`length(f_ratios) - index_max > 5`, which is the intended distance from the
most recent improvement.

## PEER

```powershell
java -jar jlinalg-0.3.5.jar confounders `
  --method peer `
  --omics expression.tsv `
  --pheno phenotype.tsv `
  --id IID `
  --full-design age,sex,ancestry1,ancestry2 `
  --factors 20 `
  --add-mean `
  --max-iterations 1000 `
  --seed 7201 `
  --out results/expression-peer
```

The dense variational updates follow official PEER 1.3: Gaussian factor and
loading posteriors, Gamma ARD precisions, feature-specific noise precisions,
known-covariate precision 100, and the original default priors
`alpha=(0.001,0.1)` and `epsilon=(0.1,10)`. The seed is explicit because the
official implementation initializes hidden factors and loadings randomly.

Inspect `converged`, the evidence-bound history available from the Java API,
and the residual-variance path. Factors can capture ancestry or protected
biology if those variables are not supplied as observed covariates.

## ComBat

```powershell
java -jar jlinalg-0.3.5.jar batch-adjust `
  --method combat `
  --omics expression.tsv `
  --pheno phenotype.tsv `
  --id IID `
  --batch sequencing_plate `
  --preserve case,age,sex `
  --out results/expression-combat
```

The default is parametric empirical Bayes. `--nonparametric`, `--mean-only`,
and `--reference-batch VALUE` expose the corresponding Bioconductor contracts.
A singleton batch forces mean-only adjustment. Features with zero variance in
any non-singleton batch are returned unchanged and counted on screen.

ComBat rejects rank-deficient covariates and covariates confounded with batch.
This is intentional: when batch and biology cannot be separated, software must
not invent an adjusted answer. Plain ComBat is not appropriate for raw counts;
normalize appropriately or use a separately validated count model.

## Use factors in association models

The default scientific workflow is to merge `*.factors.tsv` into the phenotype
table and include the factors explicitly:

```powershell
java -jar jlinalg-0.3.5.jar `
  --pheno phenotype-with-svs.tsv `
  --id IID `
  --omics expression.tsv `
  --formula "trait ~ age + sex + SV1 + SV2 + <omics>" `
  --out results/association.tsv
```

Do not estimate factors once on the complete dataset and reuse them across
prediction-validation folds. Frozen projection and fold-owned preprocessing
remain a separate workflow because refitting on held-out samples leaks their
distribution into training.

## Reproduce the reference checks

The frozen fixture was generated from Bioconductor `sva` 3.60.0 at commit
`87a4798db8134fd0952a516489dd555545a92d36` and the AutoSVA source file:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' `
  src/test/R/confounder-reference.R `
  build/reference-sources/sva-3.60.0 `
  'D:/git/NIH-R/src/sva-lite.R' `
  src/test/resources/r-reference/confounding

.\gradlew.bat test `
  --tests org.jlinalg.confounding.ConfounderRReferenceTest
```

See [the validation report](../latent-confounder-validation.md) for numerical
criteria, source provenance, performance choices, and unimplemented boundaries.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Jeffrey T. Leek and John D. Storey (2007) — Capturing heterogeneity in gene expression studies by surrogate variable analysis](../CITATIONS.md#leek-storey-2007)
- [Oliver Stegle et al. (2012) — Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses](../CITATIONS.md#stegle-peer-2012)
- [W. Evan Johnson, Cheng Li, and Ariel Rabinovic (2007) — Adjusting batch effects in microarray expression data using empirical Bayes methods](../CITATIONS.md#johnson-combat-2007)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
