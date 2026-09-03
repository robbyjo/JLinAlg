# TOPMed 100-gene GAM performance

Profiled on 2026-09-03 using 5,027 complete observations and the first 100
eligible, nonconstant genes. The fitted model was:

    BMI ~ s(omics, P-spline k=10, second-difference penalty)
        + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred
        + EO_PER_Pred + BA_PER_Pred

R used mgcv::gam with bs="ps", k=10, m=2, and method="REML". JLinAlg used
PreparedGamPredictorScan and exact Gaussian coefficient-space REML. Data
loading, the one-time compute-context selection, and one warm-up gene are
excluded. Each value below is the median of three complete 100-gene scans.
Native math threads were fixed at one; JLinAlg parallelism is across genes.

## Results

| Runtime | Workers | Median seconds | Genes/second | Speed vs R |
|---|---:|---:|---:|---:|
| R / mgcv | 1 | 16.620 | 6.02 | 1.00x |
| JLinAlg / preferred | 1 | 1.160 | 86.20 | 14.33x |
| JLinAlg / preferred | 8 | 0.337 | 297.09 | 49.38x |

Eight workers delivered a 3.45x speedup over one JLinAlg worker. The first
parallel measurement was slower because only one worker had been JIT-warmed;
the median retains the middle complete-scan observation rather than reporting
the best run.

## Numerical checks

All 100 JLinAlg fits were finite. One-worker and eight-worker results were
bit-for-bit identical for smooth EDF, smoothing parameter, likelihood, fitted
checksum, and residual sum of squares.

JLinAlg and mgcv use slightly different P-spline knot placement and penalty
scaling, so their raw smoothing parameters and restricted-likelihood constants
are not directly comparable. Across the 100 fits:

- median absolute smooth-EDF difference: 0.0692;
- maximum absolute smooth-EDF difference: 2.0066;
- median relative RSS difference: 0.00812%;
- maximum relative RSS difference: 0.13312%.

mgcv warned that some basis coefficients contained no information for three
genes per measurement. Those genes were retained in both runtimes. A separate
unit regression compares the new coefficient-space path against JLinAlg's
former dense observation-space REML and matches fixed effects, variance
components, fitted values, and random-effect prediction-error variances.

## Optimization

The previous Gam.fitGaussian implementation materialized an n-by-n smooth
covariance and repeatedly factorized it. A single real TOPMed gene exceeded
60 seconds in that path.

The optimized implementation:

- keeps a 10-basis P-spline as eight penalized random coefficients;
- performs REML through Woodbury/Henderson coefficient-space equations;
- computes X'X, X'y, Z'Z, Z'X, Z'y, and y'y once before
  smoothing-parameter optimization;
- recovers small-system prediction-error variances so smooth EDF remains
  available;
- retains one backend context across a changing-predictor scan;
- supports bounded parallel fitting across genes while native BLAS remains
  single-threaded.

## Reproduction

    & 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' src\benchmark\r\topmed_gam_benchmark.R --prepared_dir build\benchmarks\topmed100 --genes 100 --measurements 3

    $env:MKL_NUM_THREADS='1'
    $env:OMP_NUM_THREADS='1'
    .\gradlew.bat benchmarkTopmedGam --no-daemon '-PtopmedGamArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 1 --measurements 3 --backend preferred --output-prefix build/benchmarks/topmed100/jlinalg_gam_preferred_t1'
    .\gradlew.bat benchmarkTopmedGam --no-daemon '-PtopmedGamArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 8 --measurements 3 --backend preferred --output-prefix build/benchmarks/topmed100/jlinalg_gam_preferred_t8'

Timing and per-gene result CSVs are written under
build/benchmarks/topmed100.
