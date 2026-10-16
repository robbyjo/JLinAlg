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

## Draw a plot

The `plot` command reads CSV/TSV columns and writes PNG, JPEG, BMP, GIF, TIFF,
SVG, EPS, or PDF. It supports grouped scatter, line, connected-scatter, area,
step, bubble, histogram, density, dot, box, violin, bar, pie, donut, clustered
heatmap, spider/radar, circular bar, parallel-coordinate, dendrogram, treemap,
network, chord, bundled-edge, Sankey, and two-set Venn diagrams. Palettes, line widths and
patterns, marker shapes, and point sizes are configurable from the CLI.

```powershell
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type scatter --input results.csv --x estimate --y negative_log10_p `
  --group cohort --title "Association results" --out associations.pdf
```

See the [command-line plotting guide](docs/plotting.md) for mappings, themes,
axis and legend controls, missing-value behavior, and chart-specific examples.

Visit the **[JLinAlg website](https://robbyjo.github.io/JLinAlg/)** for the
complete feature guide, **[tutorials and vignettes](https://robbyjo.github.io/JLinAlg/vignettes/)**,
and Java/CLI examples, including [GRM construction](https://robbyjo.github.io/JLinAlg/vignettes/grm.html).
The [searchable scientific citation index](https://robbyjo.github.io/JLinAlg/citations.html)
links each method family to its primary publication.

## Numerical accuracy and performance

Source builds add [replicated single-cell RNA](docs/vignettes/single-cell.md) and
[physical spatial RNA workflows](docs/vignettes/spatial-analysis.md), with sparse
QC, pseudobulk, donor-aware inference, tissue graphs and spatial tests. See the
[validation and remaining scope](docs/cell-spatial-validation.md).

The source build adds [project-over-local YAML configuration](docs/vignettes/project-configuration.md),
[variant annotation, consequence adapters and evidence scoring](docs/vignettes/variant-followup.md),
and [six network workflows](docs/vignettes/network-followup.md), including WGCNA,
GENIE3 and CARNIVAL adapters. These use explicit local/remote data sources and
audited outputs. The study pipeline also offers separately opted-in AI assistance.

The source build includes the [auditable omics research pipeline](docs/vignettes/auditable-omics-pipeline.md):
a minimal versioned study contract, deterministic planning and approval,
matrix- or result-first execution, immutable amendments, hash-chained evidence,
`SESSION_HISTORY.md`, a one-command `study auto` runner, automatic public
GO/Reactome and GENCODE resource resolution through an embedded HTTPS client, a
conditional checksum-verified 1000 Genomes LD download, quantitative
cross-cohort anomaly flags, and an LLM-oriented evidence index with bounded key
findings. Repeating a command resumes verified work and retries unfinished steps;
HTML evidence and manuscript supplement previews have explicit nesting, row and
text budgets; full machine-readable evidence and checksum-linked result files
remain unchanged. Index search covers operation identity, state and summary.
An offline manuscript revision can explicitly declare `editorial_payload: {mode:
omit, rationale: ...}` to retain a full local paper without generating an external
editorial handoff. Omission is audited and prevents editorial submission; the
default 120,000-character payload limit and exact-hash approval remain unchanged.
`--rerun` starts a fresh immutable attempt. Optional data dictionaries, preapproved
sensitivity analyses, audited browser handoff (the default assistance transport),
explicit OpenAI/Gemini or OpenAI-compatible API connections, cited literature review,
and explicitly enabled HTML manuscript drafts extend the same provenance chain. See the
[hands-on tutorial](docs/vignettes/study-pipeline-tutorial.md) and
[runnable two-cohort example](examples/study-pipeline/README.md). Missing
fine-mapping, colocalization, MR, or literature inputs become explicit indexed
blockers while independent evidence branches continue. Autonomous LLM orchestration
and certified submission workflows remain future work.

For hypotheses with two or more result tables, the study pipeline automatically
creates every pairwise effect scatterplot as PNG and PDF, with the Pearson
correlation shown to four decimals in the legend. Pairs sharing a cohort ID are
labeled within-cohort phenotype comparisons; pairs from different cohorts are
labeled cross-cohort comparisons. Discovery/replication operations also create
two-set Venn diagrams. With at least three named cohorts and a completed
meta-analysis, JLinAlg automatically reruns leave-one-cohort-out replication and
meta-analysis. An omission that raises the pairwise replication fraction and
lowers mean I-squared is featured as a labeled sensitivity result; the complete
all-cohort analysis remains primary unless a prespecified rationale or reviewed
amendment says otherwise.
The [bioinformatics evidence workflow](docs/vignettes/bioinformatics-evidence.md)
adds audited public-source discovery, reviewed regional fine-mapping/colocalization/MR,
protein and reference-sequence investigation, non-LLM literature search, per-feature
candidate evidence, and local evidence archives. A separate offline synthetic example
exercises the numerical chain; missing scientific inputs remain explicit questions.

The first [study questions and decisions tranche](docs/vignettes/study-questions-and-decisions.md#implemented-first-tranche)
adds deterministic prerequisite questions, checked batch answers, answer-linked
amendment proposals, branch-specific gates, and stale-answer protection. Exact
questions and answers appear in `QUESTIONS.md`, the research notebook, and the
searchable local index. Optional LLM-suggested questions now enter the same local
review inbox; interactive prompting remains future work. See the [question workflow](docs/vignettes/study-pipeline-tutorial.md#documented-scientist-questions).

The source build adds [Tobit, AFT and ordinal likelihoods](docs/vignettes/censored-ordinal.md)
and [rare-event and survey inference](docs/vignettes/sampling-models.md).
See [regression validation](docs/regression-inference-validation.md) and the
[scientific reference page](docs/CITATIONS.md) for tested scope and assumptions.

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
