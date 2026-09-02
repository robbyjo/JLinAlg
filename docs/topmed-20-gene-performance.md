# TOPMed 20-gene performance profile

Profiled on 2026-09-02 using 5,027 complete observations and the first 20
eligible, nonconstant genes. Data loading, input preparation, and reusable
null-model setup are excluded; each reported value is the median of three
warmed measurements covering all 20 gene fits.

The benchmark uses `SampleName` as the phenotype/omics key (the file uses this
capitalization), `sabreid` as the pedigree key, `Levy_Set` as `Batch`, and adds
samples absent from the supplied pedigree as unrelated founders.

## Results

| Model | R serial (s) | CHOLMOD, 1 worker (s) | CHOLMOD, 8 workers (s) | 1 worker vs R | 8 workers vs R | Worker speedup |
|---|---:|---:|---:|---:|---:|---:|
| OLS | 0.040 | 0.008259 | 0.006290 | 4.84x | 6.36x | 1.31x |
| REML, `(1\|Batch)` | 0.860 | 0.683899 | 0.211906 | 1.26x | 4.06x | 3.23x |
| REML, pedigree + `(1\|Batch)` | 17.620 | 18.490962 | 6.230680 | 0.95x | 2.83x | 2.97x |

A comparison factor above 1 means JLinAlg is faster. With one worker CHOLMOD
is 1.26x faster for batch REML and within 5% of R for pedigree REML. Eight
workers make the two mixed-model scans 4.06x and 2.83x faster than serial R.
The very short OLS scan is dominated by sub-millisecond scheduling and timer
noise, but remains faster than R in both configurations.

The corrected native bridge was checked against R over all 60 reported fits.
The maximum absolute difference was `1.92e-05` in beta and `7.30e-06` in
standard error; one-worker and eight-worker results were bit-for-bit identical.
Multi-component pedigree fits use a bounded coordinate refinement after the
primary optimizer to avoid backend-dependent convergence from roundoff-sized
objective differences. Dedicated native tests additionally cover
factorization, log determinant, numeric refactorization with an unchanged
sparsity pattern, multi-right-hand-side solves, and CHOLMOD-first default
selection.

### JDistlib 0.10.1 CUDA comparison

The strict GPU policy selected CUDA 13.3 on an NVIDIA GeForce RTX 2080. It did
not silently route back to the CPU.

| Model | CUDA, 1 worker (s) | CUDA, 8 workers (s) | CHOLMOD/CUDA at 1 worker |
|---|---:|---:|---:|
| OLS | 0.012695 | 0.015802 | 1.54x |
| REML, `(1\|Batch)` | 2.387113 | 2.269844 | 3.49x |
| REML, pedigree + `(1\|Batch)` | >120 s setup cutoff | not run | >6.49x before any measured genes |

CUDA sparse factorization is synchronized through one backend context in this
release, so eight host workers do not materially improve it. A thread dump
during the pedigree cutoff showed the main thread waiting for completion of
the CUDA sparse Cholesky kernel during an objective refactor. This matrix is
therefore much better suited to CHOLMOD on the tested hardware. The strict
OpenBLAS policy was also probed but is unavailable on this Windows host.

## Runtime and methods

- CPU: Intel Core i9-9900K, 8 physical cores / 16 logical processors.
- JLinAlg: Java 21 and JDistlib 0.10.1; the default sparse path is
  SuiteSparse 7.12.2 /
  CHOLMOD 5.3.4 linked to oneMKL 2023.1 for BLAS/LAPACK.
- Native BLAS threads were fixed at one (`MKL_NUM_THREADS=1`,
  `OMP_NUM_THREADS=1`); JLinAlg parallelism is across genes.
- R: 4.6.1, `lme4` 2.0.6, `pedigreemm` 0.3.5, `Matrix` 1.7.5, and
  `data.table` 1.18.2.1, with all native threading variables fixed at one.
- R used `lm`, `lmer(..., REML=TRUE)`, and `pedigreemm`; JLinAlg used its fast
  OLS scan and sparse REML implementations.

The original JLinAlg comparison forced `BackendPolicy.CPU`, which is the pure
Java CPU backend: it did not use OpenBLAS, oneMKL, or the GPU. JLinAlg now uses
`BackendPolicy.PREFERRED` by default, which selects CHOLMOD when its native
bridge is present and otherwise retains the GPU/oneMKL/OpenBLAS/CPU fallback
chain. Every backend remains explicitly selectable.

## Implemented optimizations

- Analytically profile the residual variance, reducing batch REML from two
  numerical parameters to one and pedigree-plus-batch REML from three to two.
- Skip prediction covariance, random-effect recovery, fitted values, and
  residual materialization during likelihood-only optimizer evaluations.
- Retain sparse design structure, cross-product pattern, precision determinants,
  backend context, and a deterministic null-model variance start across genes.
- Retain the OLS null projection and backend context across predictor blocks.
- Use a Cholesky normal-equation fast path for full-rank OLS covariates, with the
  existing rank-aware least-squares path as fallback.
- Share one native backend context per scan, avoiding repeated initialization
  and native-context close races during parallel execution.
- Use SuiteSparse/CHOLMOD supernodal or simplicial LL' factorization with AMD
  ordering, retaining one symbolic analysis and native factor workspace per
  gene worker.

## Reproduction

Prepared inputs and result CSVs are under `build/benchmarks/topmed100`. The R
harness is `src/benchmark/r/topmed_100gene_benchmark.R`; the Java harness is
`src/benchmark/java/org/jlinalg/benchmark/Topmed100GeneBenchmark.java`.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_100gene_benchmark.R `
  --mode benchmark --prepared_dir build/benchmarks/topmed100 `
  --genes 20 --measurements 3 --models ols,reml,pedigree

$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
.\gradlew.bat nativeCholmodJar
.\gradlew.bat benchmarkTopmed100 --no-daemon `
  '-PtopmedArgs=--prepared-dir build/benchmarks/topmed100 --genes 20 --threads 8 --measurements 3 --models ols,reml,pedigree --backend cholmod --output-prefix build/benchmarks/topmed100/cholmod_t8_20_m3'

.\gradlew.bat benchmarkTopmed100 --no-daemon `
  '-PtopmedArgs=--prepared-dir build/benchmarks/topmed100 --genes 20 --threads 1 --measurements 3 --models ols,reml --backend gpu --output-prefix build/benchmarks/topmed100/gpu_t1_20_m3'
```
