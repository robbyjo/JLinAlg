# TOPMed Cox/coxme performance profile

Profiled on 2026-09-03 using the prepared TOPMed transcriptome cohort and the
first 20 genes. The fixed-effects model was:

`Surv(time, event) ~ gene + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred`

The grouped model added `(1|Levy_Set)`. The pedigree model added both a
pedigree frailty and `(1|Levy_Set)`. The generated full cohort contains 5,027
observations and 3,327 events. R 4.6.1 used survival 3.8-6 and coxme 2.2-22;
JLinAlg selected `cholmod+onemkl`. Native factorization threads were held at
one, and JLinAlg worker counts refer to independent gene fits.

## Optimized results

Times are medians of three warmed 20-gene scans. The R pedigree entry is an
explicitly labelled extrapolation from a fresh 31.32-second one-gene run; a
20-gene R pedigree scan was not repeated because it would take about ten
minutes.

| Model | R, 1 worker | JLinAlg, 1 worker | JLinAlg, 8 workers | JLinAlg vs R, 1 worker | JLinAlg vs R, 8 workers |
|---|---:|---:|---:|---:|---:|
| Fixed Cox | 0.450 s | 0.208 s | 0.052 s | 2.17x faster | 8.67x faster |
| Cox + `(1\|Levy_Set)` | 3.500 s | 1.947 s | 0.590 s | 1.80x faster | 5.93x faster |
| Cox + pedigree + `Levy_Set` | ~626.4 s (extrapolated) | 36.767 s | 10.787 s | ~17.0x faster | ~58.1x faster |

Eight JLinAlg workers produced 4.00x, 3.30x, and 3.41x speedups over the
single-worker JLinAlg scans, respectively.

For the fixed scan, the maximum absolute R/JLinAlg difference across 20 genes
was `2.43e-6` for beta and `5.62e-7` for its standard error. The single large
fixed-model difference comes from a nearly separated gene whose beta is about
-14.64; the other estimates agree more tightly. For the grouped scan, maximum
differences were `1.91e-5` for beta, `1.02e-6` for standard error, and
`2.60e-4` for the fitted frailty variance. Every fit converged.

The 100-row pedigree cross-check retained 19 nonconstant genes. JLinAlg
completed it in 0.123 seconds versus 1.220 seconds for R (9.91x faster), with
maximum beta and standard-error differences of `3.31e-4` and `1.70e-3`.

On the full-cohort one-gene pedigree fit, coxme estimated variances
`0.0101485, 0.0622674`, beta `-0.0811971`, and standard error `0.210538`.
Holding JLinAlg at those exact variances gave beta `-0.0832317` and standard
error `0.210777` in 0.093 seconds. The small remaining coefficient difference
is consistent with coxme's sparse frailty approximation rather than variance
search error; it is less than 0.01 standard errors.

## What changed

- `CoxRegression.Prepared` and `CoxMixedModel.Prepared` now own one backend
  and a cached risk-set plan for an entire scan.
- The right-censored likelihood reuses event ordering and work buffers, and
  small Newton systems use an in-process Cholesky factorization to avoid JNI
  overhead.
- Grouped Cox scans warm-start variance estimates across genes and use a
  bracketed local profile search with explicit active-bound convergence.
- `SparseCoxMixedModel` represents the pedigree through its precision matrix,
  uses grouped dense frailty columns, and solves the block system by a Schur
  complement with a reusable per-worker CHOLMOD symbolic factor.
- `PedigreeRandomEffectTerm.ofSparse` constructs Henderson's inverse
  relationship matrix directly. `ofUninbred` avoids ever materializing the
  dense 10,103-by-10,103 relationship matrix used by this data set.
- The benchmark owns a persistent worker pool and one prepared model per
  worker, so neither executors nor native factors are recreated per gene.

The sparse pedigree kernel currently requires one stratum, right-censored
data, distinct event times, and one unit-valued pedigree incidence per row.
Unsupported data should continue to use the general dense `CoxMixedModel`.

## Synthetic survival phenotype

The source phenotype has prevalent CHD/CVD indicators but no defensible
incident-event time. The R preparation script therefore generates a
deterministic right-censored phenotype with seed 20260903. Event and censor
times are independent exponentials; event hazard depends on Sex, standardized
Age and BMI, and a modest `Levy_Set` effect. Genes are absent from the
data-generating model. Continuous times make ties effectively absent, while
both runtimes use Efron handling.

## Reproduction

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_cox_benchmark.R `
  --prepared_dir build\benchmarks\topmed-obesity `
  --genes 20 --max_rows 5027 --measurements 3 --models cox,coxme

.\gradlew.bat benchmarkTopmedCox --no-daemon `
  '-PtopmedCoxArgs=--prepared-dir build/benchmarks/topmed-obesity --genes 20 --max-rows 5027 --threads 1 --measurements 3 --models cox,coxme,pedigree --backend preferred'

.\gradlew.bat benchmarkTopmedCox --no-daemon `
  '-PtopmedCoxArgs=--prepared-dir build/benchmarks/topmed-obesity --genes 20 --max-rows 5027 --threads 8 --measurements 3 --models cox,coxme,pedigree --backend preferred'
```
