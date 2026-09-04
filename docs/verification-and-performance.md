# Accuracy and performance verification

JLinAlg 0.2.0 checks every shipped component family with automated correctness
tests and an executable performance route. Independent cross-language
fixtures are used whenever a directly comparable implementation exists.
“Performance checked” means that a deterministic benchmark or representative
profile was run; it does not mean that every API is universally faster than
every alternative.

The clean v0.2.0 release gate discovered 253 tests: 250 passed, zero failed,
and three optional CHOLMOD-native tests were skipped because the native
library was not staged in that portable build. Separate native benchmark and
parity runs are documented for CHOLMOD.

## Accuracy matrix

| Component family | Independent or contract reference | What is checked |
| --- | --- | --- |
| OLS and GLM | base R `lm`/`glm`, statsmodels | coefficients, covariance, tests, likelihood/deviance |
| REML and mixed models | `nlme`, `lme4`, `rrBLUP`, `pedigreemm` | fixed effects, variance components, BLUP/PEV, sparse/dense parity |
| GEE | `geepack` and `geer` | Gaussian/binomial/Poisson coefficients, robust covariance, working correlation |
| GAM/GAMM | `mgcv` and `lme4` | fitted values, EDF, likelihood, smooth and random-effect quantities |
| Distributional/vector models | `gamlss` and `VGAM` | parameter predictors, likelihoods, probabilities |
| Cox and frailty | R `survival` and `coxme` | coefficients, SEs, likelihood, grouped and pedigree frailty |
| Penalized regression | `glmnet` plus direct normal equations | paths, objectives, deterministic CV, ridge solution |
| Association and set tests | scalar/batch parity and GMMAT | OLS/GLM/REML scans, Burden, SKAT, SKAT-O |
| Mendelian randomization | `MendelianRandomization` and base-R equations | IVW/Egger/median, heterogeneity, LD, xWAS screening |
| Meta-analysis | base-R matrix/profile calculations | fixed/random estimates, heterogeneity, meta-regression |
| Time series | R `stats::arima` and deterministic identities | ARIMA likelihood, forecasts, selection, diagnostics |
| SuSiE and colocalization | `susieR` and `coloc` | PIP, coefficients, ELBO, credible sets, H0-H4 posteriors |
| Structural equation models | `lavaan` 0.7-2 | paths, variances/covariances, SEs, likelihood and fit indices |
| Compute backends and file pipelines | cross-backend parity and format fixtures | dense/sparse solves, genotype/omics parsing, deterministic parallel output |

The full Java verification suite runs without R because expected values and
source versions are committed. Regeneration scripts live under
`src/test/resources/r-reference` and `src/benchmark/r`.

## Published measured results

The table reports real medians from the linked workload documents on the
project's Windows 11 Intel Core i9-9900K development host. “JLinAlg” is the
fastest documented production configuration for that row, commonly eight
outer workers; native/R threads and timed-region exclusions are specified in
each source document.

| Workload | Reference | Reference time | JLinAlg time | Measured speedup |
| --- | --- | ---: | ---: | ---: |
| 20-gene OLS scan | R `lm` | 0.040 s | 0.006290 s | 6.36x |
| 20-gene batch REML | R `lmer` | 0.860 s | 0.211906 s | 4.06x |
| 20-gene pedigree REML | R `pedigreemm` | 17.620 s | 6.230680 s | 2.83x |
| 20-gene grouped Laplace GLMM | R `glmer` | 18.160 s | 0.556 s | 32.69x |
| Three 100-lambda penalized paths | R `glmnet` | 0.086000 s | 0.020893 s | 4.12x |
| 100-gene Gaussian GAM scan | R `mgcv` | 16.620 s | 0.337 s | 49.38x |
| 100-gene batch GAMM scan | R `mgcv` | 21.940 s | 0.559 s | 39.24x |
| 20-gene exchangeable GEE | R `geer` | 562.530 s | 0.169 s | 3332.2x |
| 20-gene fixed Cox scan | R `coxph` | 0.450 s | 0.052 s | 8.67x |
| 223,963 random-effects meta-analyses | vectorized base R | 8.470 s | 0.159 s | 53.2x |
| 45,000-pair xWAS MR screen | base-R loop | 0.310000 s | 0.052365 s | 5.92x |
| Ten-window BMI SKAT-O suite | R `GMMAT` | 42.220 s | 7.899 s | 5.35x |
| SuSiE, 574 by 1,001 | `susieR` | 0.340 s | 0.0704391 s | 4.83x |
| SEM, 4,680 rows and 36 parameters | `lavaan` | 0.0700000 s | 0.0081623 s | 8.58x |

These ratios are not interchangeable: some time complete scans, others one
fit, and each excludes only the preparation described in its source. The
unusually large GEE ratio reflects large clusters and JLinAlg's structured
equicorrelation solve rather than a general 3332x promise.

## Reproduce and interpret

Exact commands, package versions, model contracts, numerical differences, and
limitations are recorded in:

- [TOPMed 20-gene OLS/REML](topmed-20-gene-performance.md)
- [TOPMed obesity GLMM](topmed-obesity-glmm-performance.md)
- [TOPMed penalized regression](topmed-penalized-performance.md)
- [TOPMed GAM](topmed-gam-performance.md) and
  [GAMM](topmed-gamm-performance.md)
- [TOPMed GEE](topmed-gee-performance.md)
- [TOPMed Cox](topmed-cox-performance.md)
- [TOPMed meta-analysis](topmed-meta-analysis-performance.md)
- [xWAS MR](xwas-mr-cli-performance.md)
- [TOPMed set tests](topmed-sliding-set-performance.md)
- [SuSiE benchmarks](performance-benchmarks.md#susie-benchmark)
- [TOPMed SEM](topmed-sem-performance.md)

Run the release gate with:

```powershell
.\gradlew.bat clean check executableJar
```

Benchmark results remain environment-specific. Check numerical agreement,
convergence, estimator equivalence, and workload shape before interpreting a
speed ratio.
