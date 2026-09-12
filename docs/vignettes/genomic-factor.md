# Shared genetic factors and common-factor GWAS

Source-build feature. Fit a shared genetic factor from GWAS-derived genetic
covariances, then optionally test SNP effects on that factor. The estimator
uses the full sampling covariance of the genetic moments; ordinary SEM with
a fabricated sample size is not an equivalent analysis.

## Inputs and identification

Supply a labeled genetic covariance matrix S and a labeled sampling covariance
matrix V. The [LDSC command](ldsc.md) produces these directly. For traits A,B,C,
V rows and columns are `A:A,A:B,A:C,B:B,B:C,C:C`. Matrix labels may be reordered
in files; alignment is by label. Trait units must match across inputs.

The measurement model is `S = lambda lambda' + diag(psi)`. Factor variance is
fixed at one, the first loading is oriented positive, and residual variances
psi must be positive. The fitter minimizes
`(vech(S)-vech(S_model))' V^-1 (vech(S)-vech(S_model))`, with upper-triangle
row ordering for `vech`. It uses analytic derivatives, a damped Gauss–Newton
step and a convergence certificate based on the weighted score.

Three to sixteen traits are supported. A three-trait factor is just identified
and has zero fit-test degrees of freedom; its fit p-value is NaN. With more
traits the reported fit statistic has `t*(t+1)/2 - 2*t` degrees of freedom
under the regular asymptotic moment model. Parameter covariance is
`(J' V^-1 J)^-1`, transformed to the reported residual-variance scale.

## Fit the synthetic example

```powershell
java -jar build/cli/jlinalg-0.3.5.jar genomic-factor --s examples/xwas/factor-S.tsv --v examples/xwas/factor-V.tsv --traits A,B,C,D --out build/xwas/factor.tsv
```

`OUT` contains loadings, residual genetic variances and their standard errors.
`OUT.fit.tsv` contains fit statistic, degrees of freedom, p-value and iterations.
`OUT.parameter-covariance.tsv` contains the full reported parameter covariance.
`OUT.metadata.tsv` records identification and inference boundaries.

To use LDSC output instead:

```powershell
java -jar build/cli/jlinalg-0.3.5.jar genomic-factor --s build/xwas/ldsc.tsv.S.tsv --v build/xwas/ldsc.tsv.V.tsv --traits A,B,C,D --out build/xwas/ldsc-factor.tsv
```

S must be positive semidefinite and V positive definite. Negative genetic
variances, singular sampling information, nonconvergence and near-boundary
residual variances reject with no published estimates. No nearest-PD or ridge
repair is performed. A failed fit can indicate insufficient genetic signal,
incompatible traits, too few independent blocks, or a poor factor model.

## Add a common-factor SNP scan

The GWAS table has `variant ea oa beta_A se_A beta_B se_B ...`. All beta
columns must refer to the same row effect allele and have phenotype scales
compatible with S. Pass a labeled trait sampling-error correlation matrix R_e.
For each SNP, the command constructs `C_ij = se_i * R_e_ij * se_j`.
Use identity only when the sampling errors are independent; shared participants
usually require off-diagonal terms. This command assumes that the supplied
error correlation applies across the scanned SNPs.

```powershell
java -jar build/cli/jlinalg-0.3.5.jar genomic-factor --s examples/xwas/factor-S.tsv --v examples/xwas/factor-V.tsv --traits A,B,C,D --gwas examples/xwas/factor-gwas.tsv --sampling-correlation examples/xwas/sampling-correlation.tsv --out build/xwas/factor-scan.tsv
```

`OUT.gwas.tsv` reports unadjusted p-values and the GLS common-factor effect
`b = (lambda' C^-1 beta)/(lambda' C^-1 lambda)`, its standard error
`1/sqrt(lambda' C^-1 lambda)`, Z/p, and heterogeneity
`Q_SNP = (beta-lambda*b)' C^-1 (beta-lambda*b)` with t-1 degrees of freedom.
A large heterogeneity statistic indicates trait effects poorly described by
one common-factor effect; inspect those traits rather than interpreting only
the factor p-value.

**SNP inference is conditional on the fitted loadings.** It does not propagate
measurement-model estimation uncertainty or SNP–measurement cross-covariance.
An independently estimated, precise measurement model supports this conditional
analysis; use sensitivity analyses where loadings are uncertain. The fitted
loading covariance is provided for inspection. This is not full GenomicSEM
`userGWAS` parity, multiple-factor SEM, robust DWLS, or automatic overlap
estimation. LDSC intercepts are supplied on the z-product scale; interpreting
them as error correlations requires appropriate diagonal normalization and
assumptions about confounding and study overlap.

## Java API

```java
GenomicFactor.Fit model = GenomicFactor.fit(geneticCovariance, traits, samplingCovariance);
GenomicFactor.Association snp = GenomicFactor.associate(
    alignedTraitEffects, snpSamplingCovariance, model.loadings());
```

Import `org.jlinalg.xwas.GenomicFactor`. Matrices are row-major. The Java SNP
API also accepts a different complete C for each SNP. The CLI reads tables in
memory; use separate loci or bounded files when necessary.

Reference: [Genomic SEM paper](https://www.nature.com/articles/s41562-019-0566-x)
and [reference implementation](https://github.com/GenomicSEM/GenomicSEM).
See [independent R validation](../xwas-followup-validation.md).
