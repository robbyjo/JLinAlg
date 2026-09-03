# TOPMed penalized-regression performance

Profiled on 2026-09-03 using 5,027 complete observations, 100 eligible
nonconstant gene-expression features, and seven unpenalized covariates. The
model was:

`BMI ~ <100 genes> + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred`

The existing TOPMed preparation harness matches phenotype and omics on
case-insensitive `SAMPLENAME`, annotates genes from GENCODE v48, and matches
the pedigree on `sabreid`. It adds 515 samples absent from the pedigree as
unrelated founders. The stated penalized fixed-effect model has no pedigree
random-effect term, so the pedigree is audited during cohort preparation but
does not enter these fits.

## Results

Each model fits a 100-value descending lambda path with five warmups. A timing
sample is the average of ten complete paths, and the table reports the median
of seven samples. Centering, scaling, and covariance preparation are reused
across alpha values and excluded from the path timings for both runtimes.

| Model | R `glmnet`, 1 thread (s) | JLinAlg, 1 thread (s) | JLinAlg, 8 threads (s) | 1 thread vs R | 8 threads vs R |
|---|---:|---:|---:|---:|---:|
| Ridge | 0.023000 | 0.006635 | 0.005463 | 3.47x | 4.21x |
| LASSO | 0.037000 | 0.015643 | 0.012610 | 2.37x | 2.93x |
| Elastic net (`alpha=0.5`) | 0.045000 | 0.014530 | 0.013153 | 3.10x | 3.42x |
| All three paths | 0.086000 | 0.035018 | 0.020893 | 2.46x | 4.12x |

Single-thread JLinAlg is faster for every model and 2.46x faster for the full
three-path workload. The multithreaded suite runs independent alpha paths over
one immutable prepared design; it is 4.12x faster than serial R and 1.68x
faster than JLinAlg's serial suite. Individual path rows do not claim solver
parallelism: their small difference at eight threads is ordinary timing noise.

One-time JLinAlg preparation took 0.073 seconds in the recorded runs. It is
amortized when comparing alpha values, lambda sequences, or repeated fits.
Input CSV parsing and cohort construction are excluded everywhere.

## Numerical checks

One-thread and eight-thread JLinAlg coefficients are bit-for-bit identical.
At the terminal path value, the maximum absolute JLinAlg/R difference was
0.0412 for LASSO and 0.0242 for elastic net. Lambda endpoints agreed within
`3.7e-7` absolute.

Gaussian `glmnet` standardizes the response internally. That changes the
effective L2 part of a reported lambda after coefficients are returned to the
original response scale, so ridge coefficients at the same displayed lambda
do not represent JLinAlg's original-response ridge objective. JLinAlg ridge was
therefore checked against the direct standardized normal-equation solution;
at the benchmark tolerance (`1e-8`), the maximum absolute coefficient
difference was 0.0240. With tolerance `1e-12`, it fell below `8.4e-4`.

## Implemented optimizations

- Store standardized predictors column-major so coordinate scans are
  contiguous rather than strided through row-major storage.
- Use covariance updates for moderate `p <= n` problems, reducing each
  coordinate update from O(n) to O(p).
- Build the symmetric covariance matrix deterministically in parallel and
  specialize the common unit-weight case.
- Reuse prepared centering, scales, covariance, and response products across
  ridge, LASSO, and elastic-net paths.
- Materialize fitted values and residuals lazily; compute path RSS/objective
  from cached cross-products.
- Use glmnet-style scale-aware coordinate convergence and warm starts.
- Allow zero penalty factors for the seven unpenalized covariates and compute
  lambda maxima conditional on those covariates, including wide-matrix paths.

## Runtime and reproduction

- CPU: Intel Core i9-9900K, 8 physical cores / 16 logical processors.
- JLinAlg: Java 21; heap fixed at 1 GiB initial / 4 GiB maximum.
- R: 4.6.1 with `glmnet` 5.0 and `data.table` 1.18.2.1.
- R and native math thread variables were fixed at one.
- Prepared inputs and output CSVs are under `build/benchmarks/topmed100`.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_100gene_benchmark.R `
  --mode prepare --data-dir D:\Research\topmed\splicing-bmi\new `
  --prepared-dir build\benchmarks\topmed100 --genes 100

& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_penalized_benchmark.R `
  --warmups 5 --measurements 7 --repetitions 10 `
  --output-prefix build/benchmarks/topmed100/r_penalized_final

.\gradlew.bat benchmarkTopmedPenalized --no-daemon `
  '-PtopmedPenalizedArgs=--threads 1 --warmups 5 --measurements 7 --repetitions 10 --output-prefix build/benchmarks/topmed100/jlinalg_penalized_final_t1'

.\gradlew.bat benchmarkTopmedPenalized --no-daemon `
  '-PtopmedPenalizedArgs=--threads 8 --warmups 5 --measurements 7 --repetitions 10 --output-prefix build/benchmarks/topmed100/jlinalg_penalized_final_t8'
```
