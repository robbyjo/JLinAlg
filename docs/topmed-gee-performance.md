# TOPMed GEE performance

This benchmark uses the first 100 complete, nonconstant gene-level expression
features from the supplied TOPMed files. Matching is case-insensitive on
`SAMPLENAME`; the prepared data contain 5,027 complete observations. The
marginal Gaussian model is

```
BMI ~ gene + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred
    + EO_PER_Pred + BA_PER_Pred
```

and `Levy_Set` is the GEE cluster identifier. Both working independence and
exchangeable correlation use robust sandwich standard errors. R ran one thread
with R 4.6.1, geer 0.1.0, and data.table 1.18.2.1. Java used the CPU backend.

## Results

The Java values below are medians of three measurements after a one-gene
warm-up. R is one 20-gene measurement after the same style of warm-up; the
exchangeable R run took long enough that it was not repeated. Timings are
specific to the development machine and dataset.

| Working correlation | Genes | Runtime | Threads | Seconds | Relative to R |
| --- | ---: | --- | ---: | ---: | ---: |
| Independence | 20 | R/geer | 1 | 14.790 | 1.0x |
| Independence | 20 | JLinAlg | 1 | 0.222 | 66.6x faster |
| Independence | 20 | JLinAlg | 8 | 0.042 | 349.9x faster |
| Exchangeable | 20 | R/geer | 1 | 562.530 | 1.0x |
| Exchangeable | 20 | JLinAlg | 1 | 0.655 | 859.3x faster |
| Exchangeable | 20 | JLinAlg | 8 | 0.169 | 3332.2x faster |
| Independence | 100 | JLinAlg | 1 | 0.930 | n/a |
| Independence | 100 | JLinAlg | 8 | 0.177 | n/a |
| Exchangeable | 100 | JLinAlg | 1 | 3.205 | n/a |
| Exchangeable | 100 | JLinAlg | 8 | 0.819 | n/a |

The eight-thread speedup over JLinAlg single-thread was 5.26x for independence
and 3.88x for exchangeable on the 20-gene comparison. All 100 Java fits
converged. R/geer completed every independence fit but rejected two of the 20
exchangeable fits because their estimated correlation fell just outside the
positive-definite lower bound. JLinAlg constrains exchangeable correlation to
its admissible range and completed those fits.

For the 20 genes, all 20 independence fits agreed with R to maximum absolute
differences of 4.97e-14 in the gene coefficient and 1.47e-11 in its robust
standard error. Across the 18 exchangeable R fits that completed, the maxima
were 8.15e-7 and 2.58e-7, respectively.

The speed difference is driven by the workload shape and the optimized scan
path. `Levy_Set` has only four clusters, with 1,093 to 1,504 observations per
cluster. JLinAlg uses diagonal or equicorrelation identities for linear-time
working-covariance solves and computes exchangeable pair products from sums;
it does not construct or factor a dense 1,504 by 1,504 matrix. The prepared
scan also reuses cluster ordering, the backend context, and a null-model start,
parallelizes across genes, and omits full-fit residual, QIC, and influence
objects that a coefficient scan does not consume.

## Statistical limitation

This is a useful large-cluster performance stress test, but four clusters are
not enough for dependable cluster-robust inference with nine mean-model
coefficients. In particular, the sandwich meat has rank at most four, and
cluster-t inference requiring `clusters - parameters` degrees of freedom is
undefined. The reported coefficient and standard-error agreement validates
implementation parity, not scientific validity. A production analysis should
use a defensible identifier yielding substantially more independent clusters,
or use a model whose dependence structure matches the study design.

## Reproduction

Prepare the shared inputs once:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/benchmark/r/topmed_100gene_benchmark.R `
  --mode prepare `
  --data_dir D:/Research/topmed/splicing-bmi/new `
  --prepared_dir build/benchmarks/topmed100 `
  --genes 100
```

Run the 20-gene R comparisons separately so each has its own outputs:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/benchmark/r/topmed_100gene_benchmark.R `
  --mode benchmark --prepared_dir build/benchmarks/topmed100 `
  --genes 20 --measurements 1 --models gee-independence `
  --output_prefix r_independence_20

& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/benchmark/r/topmed_100gene_benchmark.R `
  --mode benchmark --prepared_dir build/benchmarks/topmed100 `
  --genes 20 --measurements 1 --models gee-exchangeable `
  --output_prefix r_exchangeable_20
```

Run Java with one or eight outer predictor workers, changing `--genes` to 20
for direct comparison with R:

```powershell
.\gradlew benchmarkTopmed100 `
  '-PtopmedArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 1 --measurements 3 --models gee-independence,gee-exchangeable --backend cpu --output-prefix build/benchmarks/topmed100/jlinalg_gee_t1'

.\gradlew benchmarkTopmed100 `
  '-PtopmedArgs=--prepared-dir build/benchmarks/topmed100 --genes 100 --threads 8 --measurements 3 --models gee-independence,gee-exchangeable --backend cpu --output-prefix build/benchmarks/topmed100/jlinalg_gee_t8'
```
