# TOPMed sliding-window rare-variant set-test profile

Profiled on 2026-09-03 using 4,761 complete, genotype-matched observations
and ten overlapping 2 kb windows advanced in 500 bp steps. Only variants
with matched-cohort MAF at or below 1% were retained. Matching, pedigree
construction, annotation overlap, MAF filtering, allele orientation, missing
dosage imputation, Beta(1,25) weighting, null-model fitting, and result export
were all completed before or after the measured analysis region.

The two fixed-effect null models were:

```text
BMI     ~ Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred
Obesity ~ Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred
```

Genotypes and phenotypes were matched by `framid`. Pedigrees used `sabreid`;
matched samples absent from the pedigree file were represented as singleton
founders. Both null models included the pedigree additive relationship.

## Results

Times are medians of three warmed measurements, each covering all ten windows.
The speed ratio is R time divided by JLinAlg time, so a value above one favors
JLinAlg.

| Outcome | Analysis | JLinAlg oneMKL (s) | R GMMAT (s) | Speed ratio |
|---|---|---:|---:|---:|
| BMI | Burden | 5.003 | 40.390 | 8.07x JLinAlg |
| BMI | SKAT | 11.689 | 40.390 | 3.46x JLinAlg |
| BMI | SKAT-O | 53.980 | 41.260 | 1.31x R |
| BMI | shared suite | 60.799 | 42.220 | 1.44x R |
| Obesity | Burden | 4.971 | 43.710 | 8.79x JLinAlg |
| Obesity | SKAT | 11.529 | 44.890 | 3.89x JLinAlg |
| Obesity | SKAT-O | 54.250 | 45.890 | 1.18x R |
| Obesity | shared suite | 60.494 | 45.780 | 1.32x R |

The score projection is the dominant cost for Burden and SKAT. Materializing
the retained projection once and multiplying each window through oneMKL made
JLinAlg 8.1--8.8x faster for Burden and 3.5--3.9x faster for SKAT. SKAT-O has
a deliberately different calibration cost: JLinAlg runs a deterministic
10,000-draw correlated-null simulation, while GMMAT uses its analytic
SKAT-O calculation. R is therefore 1.18--1.31x faster for SKAT-O and
1.32--1.44x faster for the suite on this configuration. These rows compare the
shipped methods, not identical SKAT-O calibration algorithms.

## Numerical checks

Burden and SKAT closely track GMMAT. Agreement below is calculated across all
ten windows using `-log10(p)`; the difference column is the median absolute
difference on that scale.

| Outcome | Analysis | Pearson correlation | Median absolute difference |
|---|---|---:|---:|
| BMI | Burden | 1.000000 | 0.00000124 |
| BMI | SKAT | 0.999998 | 0.000425 |
| Obesity | Burden | 0.999261 | 0.0628 |
| Obesity | SKAT | 0.999750 | 0.1441 |
| BMI | SKAT-O | 0.957748 | 1.2001 |
| Obesity | SKAT-O | 0.984193 | 0.5658 |

The larger SKAT-O differences are expected from the simulation-versus-analytic
calibration, including JLinAlg's Monte Carlo p-value floor of `1 / 10001`.
The smaller binary-model shifts reflect the independently fitted first-order
Laplace/PQL-style null in JLinAlg versus GMMAT's binary mixed-model null.

## Selected windows

All selected windows are on `chr22` and overlap
`ENSG00000294541|ENSG00000294541.1`. The ten-window union contains 486 unique
rare variants.

| Window | Start | End | Rare variants |
|---|---:|---:|---:|
| window01 | 10,522,501 | 10,524,500 | 41 |
| window02 | 10,523,001 | 10,525,000 | 38 |
| window03 | 10,523,501 | 10,525,500 | 67 |
| window04 | 10,524,001 | 10,526,000 | 165 |
| window05 | 10,524,501 | 10,526,500 | 255 |
| window06 | 10,525,001 | 10,527,000 | 308 |
| window07 | 10,525,501 | 10,527,500 | 330 |
| window08 | 10,526,001 | 10,528,000 | 250 |
| window09 | 10,526,501 | 10,528,500 | 168 |
| window10 | 10,527,001 | 10,529,000 | 113 |

## Runtime and implementation

- Host: Intel Core i9-9900K, 8 cores / 16 logical processors, 64 GiB RAM.
- JLinAlg: Java 21, JDistlib 0.10.1, strict oneMKL 2023.1 backend with eight
  native threads. The installed NVIDIA RTX 2080 was not selected because the
  dense FP64 workload ran on the explicitly requested oneMKL path.
- R: 4.6.1 and GMMAT 1.5.0, with native thread variables fixed at one.
- BMI null: pedigree sparse REML in JLinAlg (genetic variance 13.29745,
  residual variance 16.38419) and Gaussian `glmmkin` in R.
- Obesity null: sparse first-order Laplace working model in JLinAlg (pedigree
  variance 0.64833) and binomial-logit `glmmkin` in R.

The JLinAlg implementation exposes reusable prepared variant sets, generic
continuous/binary score-null projections, and a shared Burden/SKAT/SKAT-O
suite. It uses a symmetric eigendecomposition for positive-semidefinite score
covariances and batches all SKAT-O null draws into an accelerated matrix
multiply. Filtering is exposed separately so callers can keep preparation out
of timed scans.

## Reproduction

Prepared inputs, timing CSVs, p-value CSVs, and cached R null models are under
`build/benchmarks/topmed-sliding-set`. The Java harness is
`src/benchmark/java/org/jlinalg/benchmark/TopmedSlidingWindowSetBenchmark.java`;
the R harness is `src/benchmark/r/topmed_sliding_set_benchmark.R`.

```powershell
$env:MKL_NUM_THREADS='8'
$env:OMP_NUM_THREADS='8'
.\gradlew.bat benchmarkTopmedSlidingSet --no-daemon

& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_sliding_set_benchmark.R
```
