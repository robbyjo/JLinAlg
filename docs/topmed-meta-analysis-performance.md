# Four-cohort TOPMed meta-analysis performance

Profiled on 2026-09-03 using the intersection of the CARDIA, FHS, JHS, and WHI
`splicing-bmi-*geneadj.csv` exports. After omitting incomplete BMI effect/SE
rows, 223,963 transcript features had four valid cohort estimates.

Both runtimes computed inverse-variance fixed effects and profile-REML random
effects with nonnegative tau-squared. The timed region includes pooled effects,
standard errors, normal p-values, confidence and prediction intervals, Cochran
Q and its p-value, tau-squared, I-squared, and H-squared. File loading, joining,
input validation/preparation, CSV output, garbage collection, and one warm-up
fit are outside the stopwatch. R 4.6.1 used a vectorized base-R matrix kernel
with `data.table` 1.18.2.1 for untimed input. Java used Temurin 21.0.4. Native
math threads and R data-table threads were fixed at one. The host was an
8-core/16-thread Intel Core i9-9900K.

## Results

Times are medians of five complete scans.

| Model | R vectorized, 1 thread | JLinAlg, 1 thread | JLinAlg, 8 threads | Java speedup, 1 thread | Java speedup, 8 threads |
|---|---:|---:|---:|---:|---:|
| Fixed effect | 0.130 s | 0.0730 s | 0.0172 s | 1.78x | 7.56x |
| Random effects, profile REML | 8.470 s | 0.826 s | 0.159 s | 10.3x | 53.2x |

The previous general scalar API was also measured on the first 1,000 real
features. Its median throughput was about 76,147 fixed-effect and 29,017 REML
fits/second. The full single-thread batch sustained about 3,067,419 and 271,293
fits/second, respectively. This comparison includes the old path's repeated
study/result allocation and backend selection and shows where the batch API's
gain comes from.

## Numerical checks

All 447,926 final Java model rows were matched to R by model and transcript.
Maximum absolute differences were `1.85e-4` for beta, `2.32e-4` for standard
error, `3.03e-7` for p-value, `5.68e-13` for Cochran Q, and `8.53e-14` for
I-squared. The maximum scale-adjusted beta and standard-error differences were
`4.23e-6` and `2.23e-6`.

Tau-squared had a maximum absolute difference of `10.01` and maximum
scale-adjusted difference of `9.53e-5`. These cases have unusually large,
nearly flat REML profiles; the much smaller pooled-effect and standard-error
differences above describe their inferential impact. On a separate 1,000-row
direct comparison with the established JLinAlg scalar implementation, maximum
scale-adjusted differences were `9.29e-8` for beta, `1.29e-7` for standard
error, `2.57e-8` for p-value, and `1.27e-6` for tau-squared.

## Optimization

`PreparedMetaAnalysisBatch` replaces repeated tiny general linear-algebra fits
with intercept-only sufficient statistics. It:

- validates, copies, and squares standard errors once;
- keeps effects and variances in contiguous row-major arrays;
- evaluates fixed, DerSimonian-Laird, Paule-Mandel, and profile-REML models
  without allocating study objects or one-column design matrices;
- preserves the existing bounded profile-REML search, including the boundary
  comparison, rather than substituting a different tau estimator;
- computes result fields into columnar arrays; and
- partitions independent analyses into deterministic chunks for bounded
  caller-selected parallelism.

The existing `MetaAnalysis.fit` API and general meta-regression backend remain
unchanged.

## Reproduction

```powershell
.\gradlew.bat benchmarkTopmedMetaAnalysis `
  '-PtopmedMetaArgs=--input-dir D:/Research/topmed/splicing-bmi/new --measurements 5 --threads 1'

.\gradlew.bat benchmarkTopmedMetaAnalysis `
  '-PtopmedMetaArgs=--input-dir D:/Research/topmed/splicing-bmi/new --measurements 5 --threads 8'

& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_meta_analysis_benchmark.R --measurements 5
```

The prepared four-cohort intersection, timing CSVs, and per-feature result CSVs
are written under `build/benchmarks/topmed-meta-analysis`.
