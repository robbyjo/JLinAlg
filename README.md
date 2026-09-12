# JLinAlg

JLinAlg is a Java library and command-line toolkit for statistical analysis,
including regression, mixed models, omics association, and meta-analysis.
It uses [JDistlib](https://github.com/robbyjo/JDistlib) for numerical computation.

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

The self-contained executable is `build/cli/jlinalg-0.3.5.jar`.
[Published releases](https://github.com/robbyjo/JLinAlg/releases) are also
available; newer source-build features may not be in the latest release asset.

## Run an analysis

Fit a linear regression using the included synthetic phenotype table:

```shell
java -jar build/cli/jlinalg-0.3.5.jar --pheno examples/quickstart/phenotype.tsv --id sample --formula "trait ~ age" --out build/quickstart/results.tsv
```

The result contains coefficient estimates, standard errors, and p-values,
with a companion run log. Choose a fresh output path for another run.

Visit the **[JLinAlg website](https://robbyjo.github.io/JLinAlg/)** for the
complete feature guide, **[tutorials and vignettes](https://robbyjo.github.io/JLinAlg/vignettes/)**,
and Java/CLI examples, including [GRM construction](https://robbyjo.github.io/JLinAlg/vignettes/grm.html).

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

The [estimator extension guide](docs/estimator-extensions.md) covers mixed and
ordinal SEM, DWLS/WLSMV, multigroup invariance, ARIMA regression and historical
smoothing, mixed-type kernels, quantile inference, and non-tensor integration.

## License

[GNU General Public License, version 2 or later](LICENSE) (`GPL-2.0-or-later`).
