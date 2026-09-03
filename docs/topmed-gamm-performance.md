# TOPMed 100-gene GAMM performance

Profiled on 2026-09-03 using 5,027 complete observations, four
`Levy_Set` batch levels, a 10,103-person pedigree, and the first 100
eligible nonconstant genes. The two fitted models were:

    batch:
      BMI ~ s(omics, P-spline k=10, second-difference penalty)
          + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred
          + EO_PER_Pred + BA_PER_Pred + (1 | Levy_Set)

    pedigree:
      batch model + additive pedigree effect

The prepared analysis calls the batch column `Batch`. Individuals absent
from the source pedigree were included as singleton founders during data
preparation. An exact `pedigreemm::inbreeding` check found zero nonzero
inbreeding coefficients (maximum 0), so the memory-safe sparse benchmark's
zero-inbreeding Henderson precision is exact for this pedigree.

R batch GAMM used `mgcv::gam`, `bs="ps"`, `k=10`,
`m=2`, `s(Batch,bs="re")`, and REML. JLinAlg used
`PreparedGammPredictorScan`, sparse coefficient-space REML, CHOLMOD, and
oneMKL. Data loading, structure preparation, and one warm-up fit were excluded.
Each result is the median of three complete scans. Native math threads were
fixed at one; JLinAlg workers parallelize across genes.

The machine was an Intel Core i9-9900K (8 cores, 16 logical processors), Java
25, R 4.6.1, mgcv 1.9-4, and pedigreemm 0.3-5.

## Batch GAMM results

| Runtime | Workers | Median seconds | Genes/second | Speed vs R |
|---|---:|---:|---:|---:|
| R / mgcv | 1 | 21.940 | 4.56 | 1.00x |
| JLinAlg / preferred | 1 | 1.624 | 61.57 | 13.51x |
| JLinAlg / preferred | 8 | 0.559 | 178.86 | 39.24x |

Eight JLinAlg workers delivered a 2.90x speedup over one worker. All 100 fits
converged. Single-worker and eight-worker result CSVs were bit-for-bit
identical across smooth EDF, smoothing parameter, likelihood, fitted checksum,
RSS, and every variance component.

JLinAlg and mgcv use slightly different P-spline knot placement, penalty
scaling, and restricted-likelihood constants. Fitted-scale agreement across
the 100 genes was:

- median absolute smooth-EDF difference: 0.0808;
- maximum absolute smooth-EDF difference: 1.5633;
- median relative RSS difference: 0.00955%;
- maximum relative RSS difference: 0.15920%.

mgcv warned that some basis coefficients contained no information for three
genes per measurement. Those genes were retained in both runtimes.

## Pedigree plus batch GAMM results

| Runtime | Workers | Median seconds | Genes/second | Worker speedup |
|---|---:|---:|---:|---:|
| JLinAlg / preferred | 1 | 86.089 | 1.16 | 1.00x |
| JLinAlg / preferred | 8 | 23.662 | 4.23 | 3.64x |

All 100 fits converged in every measurement. The single-worker and
eight-worker per-gene outputs were bit-for-bit identical for every reported
quantity.

A fair R pedigree-GAMM timing is not reported. The installed R stack has
mgcv, lme4, pedigreemm, and coxme, but not a sparse fitter that jointly
estimates the P-spline variance, batch variance, and additive pedigree
covariance. Expressing the 10,103-level pedigree as an mgcv random-effect
smooth would densify a coefficient design wider than the 5,027 observations.
Timing a fixed-spline approximation or omitting ancestors would benchmark a
different model.

## Optimization

The previous GAMM implementation built observation-scale covariance matrices.
The prepared scan instead:

- represents the smooth with eight penalized random coefficients;
- retains the 10,103-person pedigree and batch sparse structure;
- keeps one CHOLMOD symbolic factorization per worker and refreshes only
  numeric cross-products and factors for each gene;
- precomputes static precision determinants and sufficient statistics;
- solves only eight inverse columns for spline prediction-error variances,
  avoiding a full pedigree inverse while preserving smooth EDF;
- uses deterministic warm starts and a coarse likelihood prescreen for the
  changing smooth and boundary variance components;
- skips flat boundary directions, but reactivates them whenever likelihood
  probes show improvement;
- preserves bit-for-bit results across worker counts.

The prescreen also improved robustness: on the 20-gene diagnostic subset it
found higher REML likelihoods for two genes (+3.74 and +0.68) than a
single-start unconstrained BOBYQA run.

## Reproduction

    & 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' src\benchmark\r\topmed_gamm_benchmark.R --prepared-dir build\benchmarks\topmed100 --genes 100 --measurements 3

    $env:MKL_NUM_THREADS='1'
    $env:OMP_NUM_THREADS='1'
    .\gradlew.bat benchmarkTopmedGamm --no-daemon '-PtopmedGammArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 1 --measurements 3 --models batch-gamm --backend preferred --output-prefix build/benchmarks/topmed100/jlinalg_batch_gamm_t1'
    .\gradlew.bat benchmarkTopmedGamm --no-daemon '-PtopmedGammArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 8 --measurements 3 --models batch-gamm --backend preferred --output-prefix build/benchmarks/topmed100/jlinalg_batch_gamm_t8'
    .\gradlew.bat benchmarkTopmedGamm --no-daemon '-PtopmedGammArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 1 --measurements 3 --models pedigree-gamm --backend preferred --output-prefix build/benchmarks/topmed100/jlinalg_pedigree_gamm_t1'
    .\gradlew.bat benchmarkTopmedGamm --no-daemon '-PtopmedGammArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 8 --measurements 3 --models pedigree-gamm --backend preferred --output-prefix build/benchmarks/topmed100/jlinalg_pedigree_gamm_t8'

Timing and per-gene result CSVs are written under
`build/benchmarks/topmed100`.
