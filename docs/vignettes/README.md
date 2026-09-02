# JLinAlg vignettes

These worked examples start with small in-memory arrays so the statistical
model is visible. Replace those arrays with your own columnar or file-backed
data preparation; the fitting APIs consume primitive Java arrays and do not
require a dataframe runtime.

Run the complete verification suite before adapting an example:

```powershell
.\gradlew.bat check
```

## Choose a vignette

| Vignette | Features covered |
| --- | --- |
| [Linear models and GLMs](linear-models-and-glms.md) | OLS, weights, offsets, missing rows, contrasts, GLM families, ridge, LASSO, elastic net, cross-validation |
| [Additive, mixed, and distributional models](../additive-models.md) | GAM REML/GCV, multi-penalty te/ti, cyclic/thin-plate/MRF bases, GAMM PQL/Laplace, pedigree/GRM, GAMLSS/VGAM families, prepared scans |
| [REML and mixed models](reml-and-mixed-models.md) | covariance-component REML/ML, finite DF, dense/sparse LMMs, correlated effects, prediction, refit, comparison |
| [Pedigree and generalized mixed models](pedigree-and-glmm.md) | pedigree construction, dense/sparse animal models, combined random terms, GLMM PQL, pedigree GLMM PQL |
| [Cox and frailty survival models](cox-survival.md) | right-censoring, delayed entry, strata, Efron/Breslow ties, hazard ratios, Gaussian, GRM, and pedigree frailty |
| [Formulas and compute backends](formulas-and-backends.md) | model tables, contrasts, offsets, mixed formulas, nested/correlated terms, GPU/oneMKL/OpenBLAS policy |
| [Association, GWAS, and omics QTL](association-gwas-twas.md) | CSV/TSV/VCF/BCF/BGEN, cohort QC, cryptic-relatedness GRMs, fast OLS/GLM, P3D/EMMAX, omics transforms, Burden/SKAT/SKAT-O |
| [Mendelian randomization](mendelian-randomization.md) | harmonization, IVW, MR-Egger, LD, directionality, robust/outlier, multivariable, overlap and winner's curse |
| [Meta-analysis](meta-analysis.md) | fixed/random pooling, heterogeneity estimators, Knapp-Hartung, prediction intervals, meta-regression |
| [Time series](time-series.md) | AR/MA/ARMA/ARIMA/SARIMA, exact ARMA, automatic selection, diagnostics, forecasts, ARIMA-error LMM |
| [SuSiE and SEM](susie-and-sem.md) | individual/summary fine mapping, credible sets, observed-variable path models, fit indices, equality constraints |

## Common result pattern

OLS, GLM, REML/LMM, pedigree, GLMM, association, meta-analysis, and many
specialized fits expose coefficient-level effect size, SE, statistic, and
p-value. Where `associationStatistics()` is available:

```java
AssociationStatistics inference = fit.associationStatistics();
double[] beta = inference.effectSizes();
double[] se = inference.standardErrors();
double[] statistic = inference.statistics();
double[] p = inference.pValues();
double[] minusLog10P = inference.negativeLog10PValues();
```

Always check a result's `converged()` flag when the estimator is iterative.
The [numerical contract](../numerical-contract.md) states the likelihood,
degrees-of-freedom, approximation, and missing-data assumptions behind each
result. The [performance guide](../performance-benchmarks.md) explains how to
measure backend and sparse/dense choices on the target machine.
