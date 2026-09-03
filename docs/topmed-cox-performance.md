# TOPMed Cox/coxme performance profile

Profiled on 2026-09-03 using the prepared TOPMed transcriptome cohort and the
first 20 genes. The Cox fixed-effects model was:

`Surv(time, event) ~ gene + Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred`

The grouped model added `(1|Levy_Set)`. The pedigree model added both
`(1|animal_id)` with numerator-relationship covariance and `(1|Levy_Set)`.
R 4.6.1 used survival 3.8-6 and coxme 2.2-22 with Efron ties. JLinAlg used
`CoxRegression` and the first-order Laplace `CoxMixedModel` with default
controls and `BackendPolicy.PREFERRED`, which selected `cholmod+onemkl`.
Native backend threads were held at one; JLinAlg worker counts refer to
independent gene fits.

## Synthetic survival phenotype

The source phenotype has prevalent CHD/CVD indicators but no defensible
incident-event time. A deterministic right-censored phenotype is therefore
generated once by `topmed_cox_benchmark.R` with seed 20260903. Event and censor
times are independent exponentials; event hazard depends on Sex, standardized
Age and BMI, and a modest `Levy_Set` effect. Genes are absent from the
data-generating model. The resulting full cohort contains 5,027 observations
and 3,327 events (66.2%). Continuous times make ties effectively absent, while
both runtimes still use Efron handling.

## Full cohort: fixed and Levy_Set frailty

Times are medians of three warmed 20-gene scans.

| Model | R, 1 worker | JLinAlg, 1 worker | JLinAlg, 8 workers | JLinAlg/R, 1 worker | JLinAlg/R, 8 workers |
|---|---:|---:|---:|---:|---:|
| Fixed Cox | 0.450 s | 14.308 s | 6.666 s | 31.80x slower | 14.81x slower |
| Cox + `(1\|Levy_Set)` | 3.500 s | 44.618 s | 13.399 s | 12.75x slower | 3.83x slower |

Eight workers yielded 2.15x speedup for fixed Cox and 3.33x for grouped Cox.
Across all 20 genes, maximum absolute R/JLinAlg differences were `4.45e-12`
for fixed beta, `1.05e-12` for fixed SE, `7.88e-6` for grouped beta, and
`4.09e-7` for grouped SE. All full-cohort JLinAlg fits converged.

## Pedigree comparison

The current JLinAlg Cox frailty implementation materializes dense random-effect
design, precision, and information matrices. A full 10,103-member pedigree is
therefore not feasible. In contrast, coxme completed a full-cohort one-gene
pedigree fit in 26.52 seconds using its block-sparse matrix implementation.
A 500-person JLinAlg pedigree probe exceeded 120 seconds and was stopped.

For an executable cross-runtime comparison, the benchmark selects 100 evenly
spaced cohort rows. Unobserved ancestors are included while constructing the
relationship matrix and then analytically marginalized to the 100 observed
animals. This gives the same observed Gaussian covariance in both runtimes.
One of the original 20 genes is constant in this subset, leaving 19 genes.

| Model, 100 rows/19 genes | R, 1 worker | JLinAlg, 1 worker | JLinAlg, 8 workers | JLinAlg/R, 1 worker | JLinAlg/R, 8 workers |
|---|---:|---:|---:|---:|---:|
| Fixed Cox | 0.070 s | 12.653 s | 6.174 s | 180.76x slower | 88.20x slower |
| Cox + `(1\|Levy_Set)` | 0.330 s | 13.033 s | 6.228 s | 39.49x slower | 18.87x slower |
| Cox + pedigree + `Levy_Set` | 1.220 s | 29.761 s | 10.254 s | 24.39x slower | 8.40x slower |

Maximum absolute differences on this subset were `2.71e-9` beta and
`6.44e-10` SE for fixed Cox, `3.45e-4` beta and `2.72e-4` SE for grouped Cox,
and `3.31e-4` beta and `1.70e-3` SE for pedigree Cox. Both reduced mixed models
estimated variances near zero. R accepted the boundary solutions; JLinAlg's
default coordinate profiler reached its iteration limit for all 19 mixed fits,
despite the close coefficient agreement.

The small-cohort elapsed times have coarse R timer resolution and should be
interpreted as throughput diagnostics rather than microbenchmarks.

## Bottlenecks and next work

1. `CoxRegression.fit` and `CoxMixedModel.fit` select and close a backend for
   every gene. On 100 rows, fixed and grouped scans take nearly identical time,
   demonstrating that backend lifecycle dominates useful computation.
2. Risk-set ordering and event-time maps are rebuilt on every Newton
   evaluation instead of being prepared once per survival response.
3. `CoxMixedModel` materializes dense combined designs, penalized information,
   random Hessians, and pedigree precision. CHOLMOD is selected but is not used
   by this dense frailty path.
4. Boundary variance convergence needs explicit active-bound handling instead
   of consuming the general coordinate iteration budget.

R parity will require a prepared Cox scan with one owned backend, cached
risk-set structure, and a sparse/block-plus-low-rank frailty solver comparable
to coxme's `bdsmatrix` representation. JDistlib already exposes the dense and
sparse factorizations needed for the first two items. The scalable Cox pedigree
solver needs an algorithmic representation change in JLinAlg rather than a new
backend primitive.

## Reproduction

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src\benchmark\r\topmed_cox_benchmark.R `
  --prepared_dir build\benchmarks\topmed-obesity `
  --genes 20 --max_rows 5027 --measurements 3 --models cox,coxme

.\gradlew.bat benchmarkTopmedCox --no-daemon `
  '-PtopmedCoxArgs=--prepared-dir build/benchmarks/topmed-obesity --genes 20 --max-rows 5027 --threads 1 --measurements 3 --models cox,coxme --backend preferred'
```
