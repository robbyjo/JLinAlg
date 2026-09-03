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
| BMI | Burden | 5.051 | 40.390 | 8.00x JLinAlg |
| BMI | SKAT | 5.374 | 40.390 | 7.52x JLinAlg |
| BMI | SKAT-O | 7.676 | 41.260 | 5.37x JLinAlg |
| BMI | shared suite | 7.899 | 42.220 | 5.35x JLinAlg |
| Obesity | Burden | 4.964 | 43.710 | 8.81x JLinAlg |
| Obesity | SKAT | 5.345 | 44.890 | 8.40x JLinAlg |
| Obesity | SKAT-O | 7.774 | 45.890 | 5.90x JLinAlg |
| Obesity | shared suite | 8.075 | 45.780 | 5.67x JLinAlg |

The optimized scan retains one backend selected and owned by the null model;
component eigensolves and SKAT-O calibration no longer perform hardware
discovery per window. Analytic SKAT-O uses GMMAT's moment-matched
one-dimensional integration structure. The historical seeded correlated-null
simulation remains available explicitly.

Relative to the original implementation, backend reuse made SKAT 2.18x faster
for BMI and 2.16x faster for Obesity. Backend reuse plus analytic calibration
made SKAT-O 7.03x and 6.98x faster, and the shared suites 7.70x and 7.49x
faster. Burden is unchanged because it did not perform component eigensolves
or SKAT-O calibration.

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
| BMI | SKAT-O | 0.999796 | 0.00631 |
| Obesity | SKAT-O | 0.998954 | 0.0756 |

SKAT-O agreement improved substantially after adopting analytic calibration.
The remaining small differences reflect Imhof component tails in JLinAlg
versus GMMAT's preferred Davies tails. Binary-model shifts also reflect the
independently fitted first-order Laplace/PQL-style null in JLinAlg versus
GMMAT's binary mixed-model null.

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
covariances and a scan-owned backend for all window operations. Analytic
SKAT-O is the default; explicitly requested simulation batches null draws into
an accelerated matrix multiply. Filtering is exposed separately so callers
can keep preparation out of timed scans.

## Reproduction

Prepared inputs, timing CSVs, p-value CSVs, and cached R null models are under
`build/benchmarks/topmed-sliding-set`. The Java harness is
`src/benchmark/java/org/jlinalg/benchmark/TopmedSlidingWindowSetBenchmark.java`;
the R harness is `src/benchmark/r/topmed_sliding_set_benchmark.R`.

```powershell
$env:MKL_NUM_THREADS='8'
$env:OMP_NUM_THREADS='8'
.\gradlew.bat benchmarkTopmedSlidingSet --no-daemon

# Optional legacy 10,000-draw calibration:
.\gradlew.bat benchmarkTopmedSlidingSet --no-daemon `
  '-PtopmedSetArgs=--skato-calibration PARAMETRIC_SIMULATION'

& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_sliding_set_benchmark.R
```
