# xWAS MR CLI and performance

The `mr-xwas` command runs the bounded, parallel exposure-by-phenotype pipeline
from long-format files. Exposure instruments must already be clumped; the same
instrument set is reused across every outcome phenotype.

## Command line

```console
java -jar jlinalg-<version>.jar mr-xwas \
  --exposure molecular-traits.clumped.tsv \
  --outcome phenotypes.tsv.gz \
  --output xwas-mr-results.tsv \
  --p-threshold 5e-8 \
  --threads 8 \
  --pair-block-size 256
```

Exactly one threshold representation is required:

- `--p-threshold X` retains `p <= X`;
- `--log10-p-threshold X` retains `log10(p) <= X`;
- `--negative-log10-p-threshold X` retains `-log10(p) >= X`.

The exposure file requires SNP, beta, SE, effect allele, other allele, and an
exposure identifier. The ID is detected as `id.exposure`, `Phenotype`, or
`gene`; override it with `--exposure-id-column`. The outcome file uses the same
association columns and `id.outcome` or `Phenotype`; override it with
`--outcome-id-column`. An optional `category` or `phenotype_category` column
stores families such as cardiovascular, kidney, or lung. Effect-allele
frequency is optional but improves palindromic-variant harmonization.

CSV, TSV, and gzip-compressed input and output are supported. Screening uses
multiplicative-random IVW by default; use `--screen-method fixed` for fixed
IVW. Full MR diagnostics are computed only for retained pairs. Relevant
controls are:

```text
--threads N
--pair-block-size N
--screen-method fixed|random
--bootstrap-replicates N
--confidence-level X
--seed N
--failures FILE
--overwrite
```

The result file contains screening beta, SE, statistic, p on three scales,
confidence interval, Q, heterogeneity p-value, I-squared, mean F, MR-Egger
slope and intercept, I-squared GX, weighted-median estimates, harmonization
exclusion counts, and warnings. A sibling failure table is always written so a
later diagnostic failure cannot silently remove a significant screen.

## Reproducible benchmark

The checked-in benchmark evaluates 300 molecular exposures against 150
phenotypes with 10 already-clumped instruments per exposure: 45,000 IVW MR
pairs. It uses deterministic synthetic effects shared by Java and R.

Both runtimes exclude data generation and initial preparation/indexing from
timing, perform two warm-up scans, and report five measured scans. The timed
stage is the first-stage multiplicative-random IVW screen. The threshold is
zero so no pair enters the more expensive hit-only bootstrap and diagnostic
stage. This measures the dominant all-pairs xWAS kernel rather than a
workload-dependent number of follow-up hits.

### Test environment

- Date: 2026-09-04
- OS: Windows 11 Pro 10.0.26200
- CPU: Intel Core i9-9900K, 8 cores / 16 logical processors, 3.60 GHz
- Memory: 64 GiB
- Java: Oracle Java 25 LTS
- R: R 4.6.1 UCRT

### Measured results

| Runtime | Threads | Median seconds | Pairs/second | Speedup vs R |
| --- | ---: | ---: | ---: | ---: |
| R base nested reference | 1 | 0.310000 | 145,161 | 1.00x |
| JLinAlg | 1 | 0.140896 | 319,385 | 2.20x |
| JLinAlg | 8 | 0.052365 | 859,353 | 5.92x |

The eight-thread median includes one slower first measured iteration; no
measurements were discarded. Raw timing CSV files are written under
`build/benchmarks/xwas-mr` when the commands below are run.

### Numerical validation

Java and R independently calculated beta, SE, statistic, p-value, Cochran Q,
heterogeneity p-value, and I-squared for the first 25 exposure-outcome pairs.
Maximum absolute differences were:

| Quantity | Maximum absolute difference |
| --- | ---: |
| beta | 2.36e-16 |
| SE | 4.72e-16 |
| statistic | 5.55e-16 |
| p-value | 5.55e-16 |
| Cochran Q | 4.44e-16 |
| heterogeneity p-value | 4.44e-16 |
| I-squared | 0 |

These differences are floating-point roundoff. The R reference implements the
same zero-intercept inverse-outcome-variance regression and multiplicative
dispersion `max(1, Q/(J-1))` used by JLinAlg. This is also the statistical
contract of the default multiplicative-random IVW method documented by
[TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/man/mr_ivw.Rd).

### Commands

Single-thread JLinAlg:

```powershell
.\gradlew.bat benchmarkXwasMr --no-daemon `
  '-Djlinalg.benchmark.parallelism=1' `
  '-Djlinalg.benchmark.output=build/benchmarks/xwas-mr/jlinalg-t1'
```

Eight-thread JLinAlg:

```powershell
.\gradlew.bat benchmarkXwasMr --no-daemon `
  '-Djlinalg.benchmark.parallelism=8' `
  '-Djlinalg.benchmark.output=build/benchmarks/xwas-mr/jlinalg-t8'
```

R reference:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/benchmark/r/xwas_mr_benchmark.R `
  --exposures 300 --outcomes 150 --instruments 10 `
  --warmups 2 --measurements 5 `
  --output build/benchmarks/xwas-mr/r-base
```

Dimensions can be changed with the corresponding R arguments or Java system
properties `jlinalg.benchmark.exposures`, `outcomes`, `instruments`,
`parallelism`, `warmups`, and `measurements`.

## Interpretation limits

The measured speedup applies to this host, synthetic dimensions, and the
all-pairs IVW screening stage. It is not a universal guarantee. File parsing,
JVM startup, LD clumping, result serialization, network access, and the
hit-dependent full diagnostic stage are outside the timed region.

The R comparison is a transparent single-thread base-R nested implementation,
not a measured installation of the TwoSampleMR package, which was unavailable
on the benchmark host. It uses the same equations and generated arrays and is
appropriate for validating the numerical kernel. Comparisons with a
vectorized, compiled, or parallel R implementation can yield a different
speed ratio.

For pipeline design, phenotype categories, multiple testing, and follow-up
diagnostics, see the [parallel xWAS MR vignette](vignettes/xwas-mr-pipeline.md).
