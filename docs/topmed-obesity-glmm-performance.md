# TOPMed Obesity GLM/GLMM performance profile

Profiled on 2026-09-02 and 2026-09-03 with the first 20 eligible genes and
5,027 complete observations (1,507 obesity cases and 3,520 controls). The
fixed-effects model was:

`Obesity ~ gene + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred`

The grouped GLMM added `(1|Levy_Set)`. The pedigree GLMM added both
`(1|Levy_Set)` and `(1|animal_id)`, with `sabreid` as the pedigree identifier
and absent pedigree members retained as unrelated founders.

R used `glm(..., family=binomial)`, `glmer(..., family=binomial, nAGQ=1)`, and
`pedigreemm(..., family=binomial)`. JLinAlg used exact binomial `Glm` and the
first-order `SparseGlmmLaplace` implementation. The latter consumes grouped
CSR designs and pedigree `A^-1` directly, owns one backend for the scan, and
keeps one reusable numerical CHOLMOD factor per association worker. Native
linear-algebra threads were fixed at one. Times are medians of three warmed
20-gene measurements.

## Sparse Laplace result

The sparse JLinAlg GLMMs are faster than R on this workload, including with a
single association worker.

| Model | R, 1 worker | JLinAlg, 1 worker | JLinAlg, 8 workers | JLinAlg/R, 1 worker | JLinAlg/R, 8 workers |
|---|---:|---:|---:|---:|---:|
| Laplace GLMM, `(1\|Levy_Set)` | 18.160 s | 2.897 s | 0.556 s | 6.27x faster | 32.69x faster |
| Pedigree Laplace GLMM | 80.550 s | 40.593 s | 11.308 s | 1.98x faster | 7.12x faster |

Eight association workers provided a 5.21x speedup for the grouped model and a
3.59x speedup for the pedigree model over the JLinAlg single-worker runs.
Every worker owns its own reusable numeric factor while all workers share the
scan-owned backend and immutable symbolic pattern.

This replaces the original dense implementation, where one grouped gene took
18.584 s and one pedigree gene did not finish within 120 seconds. The sparse
grouped scan is therefore about 128x faster than the projected dense 20-gene
time. No projected pedigree speedup is reported because the dense probe did
not complete.

## Numerical agreement

Across all 20 genes, maximum absolute differences from R were:

| Model | Beta | Standard error |
|---|---:|---:|
| Grouped Laplace GLMM | 4.31e-4 | 3.55e-4 |
| Pedigree Laplace GLMM | 1.34e-3 | 2.15e-3 |

The one- and eight-worker result files were bit-for-bit identical for all
reported betas and standard errors. Prepared scans intentionally retain a
fixed configured variance start rather than passing a fitted variance from one
gene to the next; this makes results independent of worker scheduling.

## Backend observations

`BackendPolicy.PREFERRED` selected `cholmod+onemkl`. CHOLMOD performs symbolic
analysis once per worker, then refactors the same sparse pattern as PIRLS
weights and variance components change.

JDistlib 0.10.1's CUDA sparse path also ran the grouped model, but its one-gene
measurement was 0.228 s versus approximately 0.12--0.15 s for CHOLMOD after
warm-up. The CUDA pedigree probe was manually interrupted after failing to
finish promptly and is excluded from the timing table. For this sparse,
moderate-size problem, GPU launch/transfer and sparse-solver overhead outweigh
compute savings.

The earlier exact-GLM measurements remain 13.069 s with the default backend,
6.863 s with Java CPU on one worker, and 3.298 s with Java CPU on eight workers,
versus 0.230 s in R. Exact GLM scan-owned backend and full-rank Cholesky/QR
fallback work remains separate from this sparse GLMM change. The installed
oneMKL runtime does not export `LAPACKE_dgeqp3`, so its rank-aware QR path still
needs a capability-aware fallback.

## Implementation

The sparse PIRLS mode solve factors

```text
C = Z' W Z + blockdiag(Q_k / variance_k)
```

and eliminates random coefficients with a Schur complement for the small
fixed-effect block. The Laplace objective is

```text
sum log p(y_i | eta_i) - 1/2 b' P b
  + 1/2 log|P| - 1/2 log|C|,
```

where `P = blockdiag(Q_k / variance_k)`. Both `C` and pedigree precision remain
in coefficient space; neither `A` nor an observation-scale covariance is
materialized.

The R harness is
`src/benchmark/r/topmed_obesity_glmm_benchmark.R`; the Java harness is
`src/benchmark/java/org/jlinalg/benchmark/TopmedObesityGlmmBenchmark.java`.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_100gene_benchmark.R `
  --mode prepare --data_dir D:\Research\topmed\splicing-bmi\new `
  --prepared_dir build\benchmarks\topmed-obesity --genes 20

.\gradlew.bat benchmarkTopmedObesityGlmm --no-daemon `
  '-PtopmedArgs=--prepared-dir build/benchmarks/topmed-obesity --genes 20 --threads 1 --measurements 3 --models glmm,pedigree --backend preferred --output-prefix build/benchmarks/topmed-obesity/sparse_t1'
```
