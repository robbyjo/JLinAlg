# JLinAlg

JLinAlg is a Java library and command-line toolkit for statistical analysis,
including regression, mixed models, omics association, and meta-analysis.
It uses [JDistlib](https://github.com/robbyjo/JDistlib) for numerical computation.
Rare-variant APIs and the `rare-meta` CLI include Burden, SKAT, SKAT-O,
ACAT-V, and canonical six-component ACAT-O tests.

## Requirements

- JDK 17 or newer.
- Internet access for the first build to download verified dependencies.
- CPU operation needs no accelerator software. Optional GPU backends need
  compatible drivers and runtimes; NVIDIA CUDA uses Toolkit 12.6.
  See [runtime setup](docs/runtime-requirements.md).

## Build

Use the included Gradle wrapper; no separate Gradle installation is required.

```shell
# Linux / macOS
./gradlew check executableJar
```

```powershell
# Windows
.\gradlew.bat check executableJar
```

The self-contained executable is `build/cli/jlinalg-0.3.6.jar`.
[Published releases](https://github.com/robbyjo/JLinAlg/releases) are also
available; newer source-build features may not be in the latest release asset.

## Run an analysis

Fit a linear regression using the included synthetic phenotype table:

```shell
java -jar build/cli/jlinalg-0.3.6.jar --pheno examples/quickstart/phenotype.tsv --id sample --formula "trait ~ age" --out build/quickstart/results.tsv
```

The result contains coefficient estimates, standard errors, and p-values,
with a companion run log. Choose a fresh output path for another run.

Visit the **[JLinAlg website](https://robbyjo.github.io/JLinAlg/)** for the
complete feature guide, **[tutorials and vignettes](https://robbyjo.github.io/JLinAlg/vignettes/)**,
and Java/CLI examples, including [GRM construction](https://robbyjo.github.io/JLinAlg/vignettes/grm.html).
The [searchable scientific citation index](https://robbyjo.github.io/JLinAlg/citations.html)
links each method family to its primary publication.

## Numerical accuracy and performance

Validation includes independent R comparisons and analytic/reference fixtures.
For example, tested SEM coefficient errors were below `7.4e-8`. A synthetic
20,000-variant score meta-analysis ran about **2.46× faster than RAREMETAL**
on the measured workload. These are fixture-specific results, not general
accuracy bounds or speed guarantees.

See the [accuracy and performance tables](docs/verification-and-performance.md),
[numerical audit](docs/advanced-validation.md), and
[rare-variant benchmark details](docs/rare-meta-validation.md) for methods,
measurement conditions, and limitations.
The [ACAT vignette](docs/vignettes/acat-rare-variants.md) covers stable Cauchy
tails, ultra-rare collapsing, R reference fixtures, and the fused ACAT-O path.

The [estimator extension guide](docs/estimator-extensions.md) covers mixed and
ordinal SEM, DWLS/WLSMV, multigroup invariance, ARIMA regression and historical
smoothing, mixed-type kernels, quantile inference, and non-tensor integration.

JLinAlg v0.3.6 adds [LDSC](docs/vignettes/ldsc.md),
[genetically predicted TWAS/PWAS](docs/vignettes/predicted-omics.md),
[genetic-factor GWAS](docs/vignettes/genomic-factor.md), and
[prediction score training and evaluation](docs/vignettes/prediction-scores.md),
each with CLI commands and worked examples. See their
[validation and supported scope](docs/xwas-followup-validation.md).

JLinAlg v0.3.6 also includes [binary probit and a shared GLM/GEE
prediction-and-contrast API](docs/vignettes/predictions-and-contrasts.md),
[individual-level IV/2SLS](docs/vignettes/instrumental-variable-regression.md),
[summary-only conditional-score inference](docs/vignettes/conditional-score-conditioning.md),
and scalable exact diffuse [ARIMA smoothing and parameter-aware forecasts](docs/vignettes/time-series.md).
Each guide states the estimand, validation reference, and important inference
boundaries.

JLinAlg v0.3.6 provides [PCA, standard SVA, AutoSVA, and PEER latent
factors plus ComBat batch adjustment](docs/vignettes/latent-confounders-and-batch.md).
The dedicated `confounders` and `batch-adjust` commands align feature matrices
by sample ID and write reusable factors, adjusted matrices, and audit manifests.
The five preceding additions are CLI-accessible through `--family probit`,
`glm-predict`, `iv-regression`, `conditional-score`, and `arima-regression`.

JLinAlg v0.3.6 provides [empirical-Bayes differential analysis,
coordinate-aware EWAS regions, adaptive/hierarchical multiple testing, and
multiple-imputation inference](docs/vignettes/differential-regions-testing-imputation.md).
The `differential`, `ewas-regions`, `multiple-test`, `multiple-impute`,
and `mi-pool` commands preserve complete testing families, required genome
builds, reproducible imputation streams, and auditable method sidecars. See the
[validation report](docs/inference-workflows-validation.md) for frozen
limma/voom, edgeR, and DESeq2 comparisons and explicit parity boundaries.

[Gene-set enrichment](docs/vignettes/enrichment.md) adds analysis-specific
backgrounds and selection expressions, Fisher ORA, EWAS `gsameth`, BH/BY correction,
and explicit annotation downloads. The vignette includes a runnable synthetic
example, source access requirements, and guidance for overlapping ontologies.

## License

[GNU General Public License, version 2 or later](LICENSE) (`GPL-2.0-or-later`).
