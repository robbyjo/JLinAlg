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

The v0.3.5 executable includes the `mediation`, `susie`, and `coloc` CLI
subcommands, association preflight and variant-filter updates, and recursive
pedigree identity resolution. Examples that name `jlinalg-0.3.5.jar` are
release-compatible except the explicitly marked source-build additions,
including cohort meta-analysis, rare-variant workflows, and GRM construction.

| Vignette | Features covered |
| --- | --- |
| [LDSC genetic architecture](ldsc.md) | observed-scale heritability/genetic correlation, shared block jackknife, S/V exports (source build) |
| [Genetically predicted TWAS/PWAS](predicted-omics.md) | prediction weights, allele/LD alignment, molecular Z tests and joint tissue/model inference (source build) |
| [Shared genetic factors](genomic-factor.md) | full-WLS genetic measurement model and conditional SNP effects/heterogeneity (source build) |
| [Prediction scores](prediction-scores.md) | Gaussian ridge/elastic-net training, imported polygenic weights and independent-cohort evaluation (source build) |
| [GRM tutorial and Java vignette](../grm-cli.md) | worked dosage matrix, filtering, VCF input, Java/CLI construction, interpretation, and mixed-model reuse ([website](https://robbyjo.github.io/JLinAlg/vignettes/grm.html); source build) |
| [Progressive CLI association tutorial](cli-association-tutorial.md) | phenotype-only OLS/REML, numeric omics, common and rare-variant GWAS, pedigree IDs, GLM/GLMM, Cox, and penalized regression |
| [CLI-only mediation tutorial](cli-mediation-tutorial.md) | tabular OLS mediation, common complete-case filtering, grouped REML, pedigrees, output files, and causal interpretation limits |
| [CLI-only Mendelian-randomization tutorial](cli-mr-tutorial.md) | public or custom GWAS inputs, canonical columns, LD-panel installation, clumping, harmonization, estimation, and outputs |
| [CLI-only SuSiE tutorial](cli-susie-tutorial.md) | summary statistics plus labeled or ordered LD, fine-mapping controls, PIPs, credible sets, and diagnostics |
| [CLI-only colocalization tutorial](cli-colocalization-tutorial.md) | two SuSiE effect tables, priors, overlap trimming, H0-H4, conditional shared-variant posteriors, and sensitivity checks |
| [Command-line-only workflows](command-line.md) | phenotype-only and omics OLS/GLM, exact LMM/GLMM, pedigree/GRM matching, family-qualified and singleton IDs, genotype scans, Cox, and operational controls |
| [Basic statistical tests](basic-statistical-tests.md) | unit-safe correlation/mean/variance inference, exact/asymptotic rank tests, tables, blocked designs, and ANOVA power |
| [Nonlinear models](nonlinear-models.md) | analytic Gaussian fixed means, additive mixed effects, sparse pedigree structures |
| [Mediation](mediation.md) | deterministic Gaussian path analysis, Sobel inference, ordinary and pedigree mixed effects |
| [Beta regression](../beta-regression.md) | classical mean/precision beta ML, CLI, link options, and independent R validation |
| [Linear models and GLMs](linear-models-and-glms.md) | CLI-only phenotype/omics OLS and GLM commands; Java weights, offsets, contrasts, penalized fits, and inference |
| [Regression families](regression-families.md) | multivariate OLS, multinomial logistic, supersmoother, quantile, kernel nonparametric, and partially linear semiparametric regression |
| [Advanced extensions](advanced-extensions.md) | SEM extensions, full-covariance meta-regression, sparse unstructured mixed effects, diffuse/missing time series, generalized MR, quadrature GLMM, and selection-aware penalized inference |
| [Additive, mixed, and distributional models](../additive-models.md) | GAM REML/GCV, multi-penalty te/ti, cyclic/thin-plate/MRF bases, GAMM PQL/Laplace, pedigree/GRM, GAMLSS/VGAM families, prepared scans |
| [LOESS local polynomial regression](../loess.md) | direct tricube fitting, prior weights, robust symmetric passes, prediction, leverage, R accuracy and speed |
| [REML and mixed models](reml-and-mixed-models.md) | CLI-only grouped/GRM exact REML; Java covariance components, correlated effects, prediction, refit, and comparison |
| [Pedigree and generalized mixed models](pedigree-and-glmm.md) | CLI-only pedigree LMM/GLMM; Java pedigree construction, PQL, adaptive quadrature, and ZIP/ZINB |
| [Cox and frailty survival models](cox-survival.md) | CLI-only fixed/GRM Cox; Java streamed omics, delayed entry, strata, shared frailty, and pedigree frailty |
| [Formulas and compute backends](formulas-and-backends.md) | CLI-only formula/backend preflight; Java model tables, contrasts, offsets, and nested/correlated terms |
| [Association, GWAS, and omics QTL](association-gwas-twas.md) | CLI-only numeric omics, VCF/BCF/BGEN, and GRM scans; Java prepared scans and Burden/SKAT/SKAT-O |
| [Omics transforms from the command line](omics-transforms.md) | exact `--transform` syntax, all built-in stages, pipeline recipes, missing values, failure conditions, and trusted plugin providers |
| [End-to-end Mendelian randomization](mr-end-to-end.md) | database installation, public/custom instruments, LD clumping, analysis, diagnostics, plotting exports, bidirectional and molecular MR |
| [Parallel xWAS MR](xwas-mr-pipeline.md) | reusable clumped instruments, phenotype families, bounded parallel exposure-outcome scans, all-pairs BH/FDR, scale-safe thresholds, and two-stage diagnostics |
| [xWAS MR CLI and benchmark](../xwas-mr-cli-performance.md) | long-format CLI, output schema, reproducible R validation, timing protocol and measured speedup |
| [MR estimator reference](mendelian-randomization.md) | harmonization, IVW, MR-Egger, LD, directionality, robust/outlier, multivariable, overlap and winner's curse |
| [Meta-analysis](meta-analysis.md) | omics Java/CLI examples, fixed/random pooling, heterogeneity, meta-regression, and links to rare-variant workflows |
| [Rare-variant cohort summaries and meta-analysis](../rare-variant-meta-analysis.md) | participant score export, RAREMETAL/rvtests import, single variants, equal/weighted burden, SKAT, SKAT-O, and calibration limits (unreleased) |
| [Time series](time-series.md) | AR/MA/ARMA/ARIMA/SARIMA, exact ARMA, automatic selection, diagnostics, forecasts, ARIMA-error LMM |
| [SuSiE, colocalization, and SEM](susie-and-sem.md) | individual/summary fine mapping, credible sets, multi-signal colocalization, and a compact SEM introduction |
| [Multi-signal colocalization](colocalization.md) | alignment, priors, posterior-overlap trimming, H0-H4, shared-variant posterior, diagnostics, and MR/xWAS follow-up |
| [Structural equation modeling](sem.md) | joint observed/latent RAM, means, ordinal PML, pattern FIML, robust inference, lavaan validation, and limitations |

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
