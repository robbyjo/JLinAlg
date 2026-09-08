# JLinAlg 0.3.0: remaining-code audit and release validation

This audit follows the [advanced-method audit](advanced-validation.md).
It targets remaining numerical models, shared inference, data pipelines,
and release packaging. A code audit and passing fixtures do not prove that
every possible dataset is correct; model approximations and unsupported
inference remain explicitly documented.

## Shared statistical tests and inference

Reproduced before repair:

- Multiplying both Pearson inputs by `1e100` changed correlation from
  `0.8285714285714286` to zero; `1e-100` produced infinity. Centered,
  bounded-coordinate calculations now preserve the result, including the
  corresponding t test and Fisher interval.
- Perfect positive Kendall correlation at n=30 returned p=0 instead of
  `2/30! = 7.539975257631812e-33`. Direct symmetric tails now avoid cancellation;
  a normalized recurrence avoids factorial overflow. Default exact selection
  follows R's n<50 cutoff; explicit exact requests have a documented
  50-million-cell resource bound.
- Negative and positive zero were ranked as distinct observations. They now
  tie in Kendall, Spearman, and the shared ranking/tie helpers.
- Spearman continued AS 89 beyond R's n<1290 cutoff; it now switches to t.
- Student/Welch and variance tests overflowed or underflowed when units changed
  by `1e200` or `1e-200`. Private sample copies are now rescaled before the
  JDistlib statistic and covariance calculations; reported estimates and
  intervals remain in original units.
- Classical ANOVA rejected a constant group despite positive pooled variance;
  the equal-variance path now permits it and scales all groups together.
- A NaN coefficient with zero SE returned p=0; nonfinite coefficients now have
  unavailable inference. Linear contrasts reject nonfinite/asymmetric covariance.
- Model-dimension products are checked without 32-bit wraparound. A failed
  native CHOLMOD refactor now invalidates numerical solves/log determinants
  until a positive-definite refactor succeeds. The native recovery regression
  requires the optional bridge and is not certified by a skipped portable run.
- Bootstrap result accounting now rejects malformed/nonfinite successful
  replicates; empirical SD avoids square overflow for large coefficient units.

Independent cross-review also caught constant-maximum bootstrap rounding,
an overflowing paired-t summary, a null-mean-dominated t-test scale, and
contrast-unit overflow/underflow. The repairs preserve nine identical maximum
double draws exactly, assemble t statistics without scaling sample variance by
the null mean, and normalize hypothesis rows for the Wald solve. Full covariance
PSD validation now rejects an indefinite matrix even if a selected one-row
contrast happens to have positive variance. Singular PSD covariance remains
valid for estimable contrasts. This validation adds an O(p³) spectral check;
contrast estimates/covariance can still exceed floating-point range in original
units even when the unit-invariant test statistic is representable.

`CoreStatisticsAuditTest` uses independently generated base-R fixtures across
exact and approximate rank tests, scaling extremes, and ANOVA. One hundred
random tied datasets compare the new Kendall count against an independent
quadratic oracle. `CoreInferenceAuditTest` checks invalid inference and dimensions.

### Measured Kendall performance

Both sides compute the complete two-sided asymptotic Kendall test on identical
deterministic tied inputs, including ranking and tie correction. Java uses
O(n log n) counting; base R's correlation kernel uses pair comparisons.
These are local warmed wall-clock medians, not a general R performance claim.

| Pairs | Java seconds | R seconds | R / Java | Identical tau |
| ---: | ---: | ---: | ---: | ---: |
| 5,000 | 0.0020009 | 0.130000 | 65.0 | 0.0009037736641225347 |
| 20,000 | 0.00325804 | 2.036667 | 625.1 | 0.0001043595718463651 |

The previous Java quadratic implementation took about 24 ms at 5,000 pairs
in the same development session. Java timing includes 10 warm-ups and seven
batches of ten calls; R uses three warm-ups and five batches of three calls.
Windows 11, Java 25 (Java 17 target), R 4.6.1, portable CPU. Raw R timings are
in `src/benchmark/resources/core-audit/r-timing.csv`.
The current host reports 32 logical processors and Intel Family 6 Model 183
Stepping 1. It is not the older i9-9900K host named in historical site reports.

Reproduce:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/R/core-statistics-audit.R
./gradlew.bat benchmarkCoreStatisticsAudit
```

Algorithm-selection reference: [base R cor.test documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html).

## LOESS

Uniform weights `1e-20` previously threw a no-positive-weight error, and weights
`1e200` could overflow the leverage determinant. Global prior-weight and local
moment normalization now preserve fits, leverage, and effective DF. Response
normalization prevents local RHS overflow; distance fallback handles a finite
predictor range spanning `-1e308` to `1e308`. Nearest-window search now uses
binary search rather than expanding every neighborhood one observation at a time.

The existing weighted/robust `stats::loess` fit, prediction, trace, and robust
weight fixtures still pass (tolerances `2e-12` and `2e-10`). New adversarial
tests cover common weight units from `1e-200` to `1e308`, large responses,
and overflowing predictor differences.

At n=5,000, span=0.2, degree=2, Gaussian direct surface, seven warmed median
full fits measured 23.333 ms before and 21.104 ms after (1.11x improvement),
versus R's 80 ms (3.79x). Prepared Java fits measured 20.379 ms. The three-point
fit checksum was `0.001584182714` on every path. R timing disables expensive
full statistics and uses approximate trace; Java computes diagonal leverage,
so this is not an identical inference-work comparison. No large speedup is
claimed for the normalization repair itself.

Commands: `benchmarkLoess` with `jlinalg.benchmark.warmups=5` and
`jlinalg.benchmark.measurements=7` (properties before task), and
`Rscript src/benchmark/r/loess_benchmark.R 5000 0.2 7 5`.
Raw paired medians: `src/benchmark/resources/core-audit/loess-timing.csv`.

## Remaining models and pipelines

The audit was divided into disjoint package groups, followed by independent
cross-review of different owners' repairs. The previous advanced-method fixtures
remain part of the integrated regression suite; those models were not simply
counted as newly audited again.

| Audited group | Repaired behavior | Independent evidence |
| --- | --- | --- |
| OLS, GLM, GEE, penalized fitting | Gaussian likelihood/DF, unit-sensitive rank/variance/penalties, convergence and weighted CV | [R/analytic fixtures and allocation/timing report](../src/benchmark/resources/fitting-audit-v030/evidence.md) |
| General Laplace/PQL, pedigree, dense REML | full marginal beta optimization, observed mode curvature/covariance, ancestry cycles, covariance identification/PSD | [base R, lme4, and pedigree matrix gates](../src/benchmark/resources/remaining-mixed/REPAIR.md) |
| GAM and distributional models | offsets, stable beta/NB/Poisson densities, score-qualified convergence, beta mixed LL, NB observed covariance, truthful ZI failure | [distributional/GAM audit](distributional-gam-audit.md) |
| Cox, mediation, nonlinear | risk-set cancellation and predictor units, final Jacobian and stopping, sparse assembly, Sobel scaling | [survival, lm/lmer, nls, independent frailty/pedigree references](../src/benchmark/resources/remaining-model-audit-benchmark/audit-results.md) |
| Genetics, SuSiE, coloc, MR | invalid covariance/support, weighted posteriors, stale variance state, stable MR likelihood/scaling and LD mode selection | [package fixtures and independent base-R equations](../src/benchmark/resources/genetic-audit/AUDIT.md) |
| Association, GWAS, set tests, pipeline, CLI | sample/dosage validity, convergence/row accounting, bounded disk BH, kernel scale and tail accuracy | [independent fixtures and paired benchmark report](../src/benchmark/resources/pipeline-audit/AUDIT.md) |
| Shared inference, tests, LOESS, compute wrappers | scale/tail/rank/validation defects described above | base-R fixtures, independent counting oracle and adversarial arithmetic tests |

### Representative paired measurements

All values below are local warmed milliseconds per call. Read each linked
report for rows, repetitions, package versions, convergence/error gates, and
included work. R formula/model-frame overhead is included where the package
requires it; Java generally consumes numeric designs. These are not universal
speed comparisons or identical allocation contracts.

| Workload | Java ms | R ms | Reference/work distinction |
| --- | ---: | ---: | --- |
| Weighted-offset OLS, n=3,000 | 0.523 | 0.5 | R lm; full Java diagnostics; R slightly faster |
| Weighted-offset Poisson, n=3,000 | 2.103 | 6 | R glm |
| Fixed-correlation Gaussian GEE | 4.577 | 22.5 | independently written R GLS/sandwich, not geepack |
| Weighted LASSO, three lambdas | 0.204 | 1 | glmnet path |
| Fixed beta, n=1,500 | 4.663 | 50 | betareg |
| Beta mixed, n=600/40 groups | 244.961 | 200 | glmmTMB; R faster on this corrected fit |
| ZIP, n=450 | 82.455 | 140 | Java numerical inference; R sdreport |
| ZINB, n=450 | 45.439 | 190 | Java fit-only versus R sdreport; not equal inference work |
| Start-stop Cox, Efron, n=240 | 0.346 | 2.139 | survival::coxph |
| Ordinary mediation, n=120 | 0.268 | 0.918 | three lm fits and analytic effects |
| Mixed mediation, n=120 | 1.150 | 23.521 | three lmer fits and analytic effects |
| Nonlinear fixed exponential, n=240 | 0.311 | 0.977 | stats::nls; refreshed after cross-review |
| General rare-binomial Laplace | 16.597 | 60 | lme4, plus independent scalar-mode R likelihood |
| OLS scan, 400 rows × 128 markers | 0.751 | 5.2 | cached Java projection versus fresh R QR per marker |
| Disk BH, 20,000 rows | 78.112 | 1 | Java includes disk I/O; R memory-only; not equal work |
| Rank-two kernel tail | 0.02930 | 0.022 | angular integration; R faster on this fixture |
| Rank-three kernel tail, q=100 | 0.11466 | 3.2 | independent base-R positive series, not a package timing |
| 64-SNP coloc evidence | 0.01244 | 0.2 | base-R configuration enumeration, not fresh coloc package timing |
| 32-SNP one-effect posterior | 0.02307 | 0.006 | R closed form; Java validates LD and runs IBSS; R faster |
| 41-founder LD stationary phase | 0.003293 | 0.030 | Java counts genotypes; R starts from counts; no BED I/O |

Likelihood comparisons include the full objective, not just matching displayed
coefficients. The new rare-binomial/Poisson Laplace fixtures agree with independent
R scalar-mode optimization and tightened `glmer(nAGQ=1)` likelihoods. Corrected
beta mixed LL error is below `1e-8`, superseding the old `0.122` discrepancy and
its invalid speed baseline. Coloc fixtures retain disjoint/tiny support instead
of subtracting nearly equal evidence sums. The extreme rank-two quadratic-form
review case at q=200 now agrees with its chi-square limit at approximately
`2.09e-45`; an absolute integration tolerance previously caused a 38% relative error.

Independent cross-review additionally reproduced tiny-weight Poisson SE error,
Gaussian GEE dispersion floors, a rounded positive-logit binomial tail, and an
extreme deficient-design projection failure. Predictor-aware score, likelihood,
and working precision now preserve the finite positive-tail MLE (beta
`0.17001035679439` in the frozen example), including the coupled PQL/Laplace
paths. Exact-zero Gaussian residual variance has the correct infinite likelihood
supremum. A minimum-norm SVD that cannot reproduce the verified equilibrated
projection fails explicitly rather than returning a wrong fit. Extreme CV weight
ratios remain valid when their per-fold calculation is representable.

Nonlinear convergence no longer accepts a large response-scale rounding floor:
adding `1e15` to the response/mean now reaches the near-optimal slope instead of
falsely accepting the starting value. Quantized residual plateaus still return
nonconvergence if their score cannot be certified. Near-one mediation confidence
levels use the small tail directly rather than rounding the upper probability to
one. Reviewers reran the original counterexamples, 1,400 independent tied/start-stop
Cox derivative checks, and nine nonlinear mixed returned-state checks.

Later review closed two additional numerical contracts. A rank-three mixture
at q=100 returned a tail around `1.2e-10`, exceeding its chi-square upper bound;
positive gamma-series summation now agrees with the independent R reference
near `2.5e-23`. Truncation must meet a relative bound, with an explicit 8,192-term
limit instead of an inaccurate silent fallback. Near-one critical values use
the direct lower tail and positive rank. OLS/P3D marker inference is invariant
at units `1e-160` and `1e160` when returned quantities are representable.

Duplicate, permuted/rescaled, or linearly dependent random covariance terms
are rejected before Laplace optimization. The check uses sparse precision
solve actions and a small reproducible probe set rather than allocating
observation covariance matrices. It is a conservative numerical rank check,
not a global identification theorem; nearly dependent or probe-inconclusive
models may be rejected. A separate backtransform repair preserves representable
subnormal covariance without multiplying extreme predictor scales together.

The final genetics cross-review replaced multi-start LD phase iteration with
bounded isolation of every stationary cubic root and comparison of all candidate
likelihoods and endpoints. It resolves both a valid 41-founder clump that
previously aborted and a 26-founder case whose midpoint was a nonmaximum
stationary point. Independent R `polyroot`/observed-likelihood checks agree,
and the bounded 2,000-table stress run has no failures. The same review restored
a representable SuSiE posterior mean of `1e-300` under a tiny prior, and IVW
beta/SE of `1e200` whose intermediate scale ratio had overflowed.
The final LD phase core has a separate paired base-R timing, with relative r²
error `4.22e-15`. It is not a PLINK or file-clumping benchmark: both exclude
BED I/O, and only Java includes genotype counting and reflection dispatch.

### Retained limitations and compatibility corrections

- No global optimum or general finite-sample coverage certificate is implied.
  PQL, first-order Laplace, nonlinear mixed REML linearization, diagonal sparse
  Cox random information, and Sobel zero-cross-covariance inference remain
  approximations under the documented model assumptions.
- One reproduced competing-mode ZINB refit remains uncertified. It now returns
  nonconvergence with unavailable inference, rather than false success. The
  complete reproducible counterexample and limits are in the distributional
  audit; it contributes no MLE or speed claim.
- `MultivariableMrResult.conditionalFStatistics()` never computed conditional F.
  The deprecated accessor now fails explicitly; `marginalFStatistics()` exposes
  the statistic actually available. True conditional strength needs a richer
  cross-exposure covariance model.
- MR Huber/PRESSO-style and contamination-mixture routines are not complete
  calibrated substitutes for every corresponding R package. Their dispersion,
  grid, and inference limits are retained in the MR vignette.
- Fresh coloc/susieR/geepack package timings were unavailable on this host.
  Existing package fixtures remain checked; new base-R equation references are
  labeled as such. No optional native/GPU runtime certification is inferred
  from portable CPU tests.

The [TODO inventory](../TODO.md) and [advanced-method report](advanced-validation.md)
retain the separate open extensions (ordinal WLSMV, richer covariance structures,
multidimensional quadrature, etc.). This release does not claim full R-package parity.

## Release validation

Final local integrated gate: **609 tests discovered, 606 passed, zero failures,
three optional native CHOLMOD skips**. `check benchmarkClasses assemble` passes
on the final production/test snapshot. The compiled artifacts target Java 17;
the local runtime is Java 25. Native-library discovery warnings are not native
runtime certification. The site contains 22 guides; its internal-link gate passes.

The release workflow verifies an explicit immutable version tag, builds/tests
Java 17 artifacts, and publishes executable, library, sources, Javadoc, and
SHA-256 assets. It no longer checks out an unrelated old tag merely because
the workflow file changed on main.
Published assets distinguish the executable `jlinalg-0.3.0.jar` from the thin
`JLinAlg-0.3.0-library.jar`; names differing only by case would collide on
Windows. Sources, Javadoc and `SHA256SUMS.txt` accompany them.

Reproduce the integrated build and all new executable benchmark routes:

```powershell
./gradlew.bat check benchmarkClasses assemble --no-parallel --max-workers=1
./gradlew.bat benchmarkCoreStatisticsAudit benchmarkFittingAudit benchmarkDistributionalAudit benchmarkRemainingMixedAudit benchmarkRemainingModelsAudit benchmarkGeneticAudit benchmarkPipelineOlsAudit benchmarkPipelineBhAudit benchmarkPipelineKernelAudit --no-parallel --max-workers=1
```

The numerical suite consumes checked-in independent fixtures; R is not required
to run it. Regenerating the R references requires the package versions named in
each linked report. Benchmark routes consume fitted results and assert their
accuracy/convergence contracts; their sequential integration run is a smoke
gate, while the tables above use the separately recorded warmed paired runs.

Executable CLI smoke checks cover variable-precision beta regression,
elastic-net five-fold CV, generalized LD-aware MR with native SVG output, and a
shuffled-sample omics scan. The scan preserves source order
`signal / constant / noise`, returns `ok / failed / ok`, and records
three source features, two tested features and one failure. Penalized output
leaves inferential columns blank rather than advertising naive post-selection
p-values. The website link gate and browser checks cover the vignette index and
the new nonlinear documentation; the new Java examples were separately compiled.
All nine new Gradle benchmark tasks passed in the final sequential integration
run. The existing-output CLI guard was also rechecked: a second invocation
without overwrite failed and left the result file's SHA-256 unchanged.
