# JLinAlg 0.1.0

## Unreleased

- Added `ColocSusie` multi-signal colocalization compatible with R coloc,
  including ID alignment, posterior-overlap trimming, scalar and weighted
  priors, H0-H4 summaries, and per-variant H4 posteriors.
- Added coloc 5.2.3/susieR 0.14.2 example-data regression fixtures and a
  deterministic combination benchmark.
- Exposed per-effect log Bayes factors from `SusieResult` and removed repeated
  hot-loop allocations from SuSiE IBSS updates.

JLinAlg 0.1.0 is the first tagged release of the Java statistical modeling
and statistical-genetics library.

## Highlights

- Linear, generalized linear, mixed, additive, survival, pedigree, and
  generalized mixed-model APIs with common inference results.
- Streaming GWAS, TWAS, EWAS, and PWAS inputs; prepared high-throughput scans;
  and Burden, SKAT, and SKAT-O tests.
- Mendelian-randomization harmonization, independent and LD-aware estimators,
  sensitivity analyses, multivariable MR, and overlap-aware inference.
- Meta-analysis, meta-regression, ARIMA/time-series models, SuSiE fine mapping,
  and observed-variable structural equation models.
- A CLI installer for freely available LD panels. Installed sources are
  normalized into JLinAlg's versioned, variant-major PLINK reference layout.
- A self-contained command-line JAR for Java 17 or newer.

## Performance maturity

MR, SuSiE, and SEM are implemented and tested in 0.1.0, but have not yet
received the large-workload optimization and benchmarking applied to JLinAlg's
performance-tuned model and association paths. Their APIs should be considered
initial, and users should profile them on representative data before
high-throughput use.

## Build

Run `./gradlew check executableJar` (or `.\gradlew.bat check executableJar` on
Windows). The executable artifact is
`build/cli/jlinalg-0.1.0.jar`.

See the [project website](https://robbyjo.github.io/JLinAlg/) and
[worked vignettes](https://robbyjo.github.io/JLinAlg/vignettes/) for feature
scope, assumptions, and examples.
