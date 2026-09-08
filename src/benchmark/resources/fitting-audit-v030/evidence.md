# Fitting audit after dd45e80

Reproduction: run `org.jlinalg.benchmark.FittingAuditBenchmark` (optional row
count, default 3000) and, from the repository root,
`Rscript src/benchmark/r/fitting_audit_v030.R`. R prepends `build/r-library`.
R 4.6.1, glmnet 5.0; portable Java CPU, 2026-09-08. No packages installed.
Reference-generation needs readable existing R user-library dependencies.

## Confirmed counterexamples and repairs

- Gaussian GLM on y=[1,2,5,7,9], X=[1,0:4]: slope p was 5.26e-43 instead
  of Student t p=0.00083285047; LL was -2.45647458 instead of -2.17941053.
  Fixed t/F inference, metadata, ML dispersion and rank-based parameter count.
- OLS with slope-column units 1e-16 incorrectly rejected full rank. The first
  optimized repair also failed at n=10000 with columns exactly 1e-6 and +/-1e6;
  the final fast path checks the ratio of column scales as well as extremes.
  Full-rank equilibration back-transforms coefficients and covariance; deficient
  SVD retains the original-coordinate minimum norm. Non-estimable individual
  inference is unavailable and non-estimable contrasts are rejected.
- Gamma GLM dispersion changed from .5185185 to .000466667 when y was scaled
  by 1e-8. GEE Gamma slope changed .51074603 to .49368534 and falsely converged.
  Removed response-unit variance floors; tests verify slopes, dispersion and
  covariance under unit changes. GEE covariance jitter is relative to scale.
- Grouped binomial proportions [1,2,3]*1e-14 with each trial weight 1e14:
  mean was clamped to 1e-12, beta drifted to -100.08 and failed convergence.
  Stable logit/log1p deviance now recovers mean 2e-14, beta -31.54304412,
  and analytic Fisher information; R binomial densities independently gate LL.
- Standardized LASSO slope 2.0292893 became zero after changing predictor units
  by 1e-16. Removed absolute SD cutoff and use a stable weighted norm.
- Equal weights 1e308 yielded a converged zero fit/objective. Normalize relative
  to maximum weight before summation (fitting and CV).
- A nearly exact, large-response fit returned RSS=0 while its actual residual
  RSS was 39.92495. Detect quadratic cancellation and reconstruct residual RSS.
- CV returned a selected nonconverged fit with maxIterations=1 and incorrectly
  gave folds equal weight despite unequal held-out weight totals. Reject
  incomplete paths; aggregate weighted held-out risk and grouped uncertainty.
  Reuse the initial full-data fit for automatic CV.
- Gaussian identity now needs one least-squares solve, not repeated IRLS QR;
  serial GEE accumulation no longer retains all cluster score/bread/meat arrays.

## Warm timing and CPU allocation

Seven batches, 20 fits each, after warmup; every fit must converge and contributes
to a checksum. Java allocations use the current-thread allocation counter, with
serial execution. These are per-fit allocation traffic, not live heap size.

| Workload (n=3000) | Java ms | R ms | Java bytes/fit | Java checksum over 140 fits |
|---|---:|---:|---:|---:|
| weighted offset OLS | 0.522630 | 0.5 | 569190 | 182.00698934111634 |
| weighted offset Poisson | 2.103355 | 6 | 4050434 | 224.08494210358435 |
| fixed-correlation Gaussian GEE | 4.577415 | 22.5 | 19772265 | 182.14895764808057 |
| weighted LASSO, 3 lambdas | 0.203725 | 1 | 242553 | 546.58570141596570 |
| weighted 5-fold LASSO CV | 0.862975 | not timed | 1669383 | 40.486549538030914 |

The matching R coefficient checksums were 182.00698934111639,
224.08494210358523, 182.14895764808062 and 546.58570141596567. R's Windows timer
has about 0.5 ms resolution after batching. Full Java inference/auxiliary
diagnostics are not identical-cost output contracts to raw R fits. The fixed
GEE reference is explicit R block-matrix GLS plus sandwich, not geepack.
Baseline pre-repair GEE allocated ~20,613,041 bytes/fit versus ~19,772,265 now.
Baseline CV checksum differs because its risk estimand was wrong; do not treat
its shorter timing as a valid accuracy-matched performance baseline.

## Gates and actual limits

Final isolated suite: 94 tests passed, zero skipped/failed. The final integrated
release gate reports 609 tests, 606 passed and three optional native skips,
including all six follow-up repairs. Maximum observed
Java/R errors over the frozen reference gates: coefficient 2.054e-8
(inverse-Gaussian, iterative reference tolerance), covariance 5.407e-10,
p-value 2.292e-9, log likelihood 1.137e-13, and CV MSE 7.340e-10.

Frozen fixtures cover coefficients, full covariance, p-values and likelihood
for weighted-offset OLS/Gaussian/Poisson; Gamma and inverse-Gaussian likelihood
are independently profiled. Weighted lasso/CV use glmnet; ridge and elastic net
with penalty factors [0,2] use exact two-predictor active-sign enumeration of
the same objective, avoiding glmnet's different penalty-factor normalization.
JUnit tests assert convergence before numerical comparisons. Unit-change,
rank/estimability, tiny-probability, huge-weight, RSS-cancellation and parallel
accumulation tests are separate analytic gates, not just generated happy data.

No geepack installation was available; not all GEE association/bias-correction
variants have independent R parity here. Strict backend selection is preserved,
but optional CHOLMOD skips do not certify its native runtime. Deficient SVD
fits identify rank in equilibrated coordinates and validate original-coordinate
SVD projection/leverage against that reference. Very ill-conditioned or
nonconvergent (including some wide) deficient systems are explicitly rejected;
no different identifying norm is silently substituted. Convergence checks do
not establish identification or rule out separation. Arbitrary floating-point
overflow and unrepresentable probability tails remain limits. Custom-family
likelihoods retain their supplied-dispersion contract. Gamma AIC intentionally
uses its profiled density, not stats::glm's approximate Gamma AIC expression.

## Independent cross-review follow-up: six repaired gates

The independent worker reproduced and rechecked these cases in
`build/meitner-cross-review/ReviewProbe.java` (raw repaired output in the same
directory). Persistent regressions are the owned `*EdgeAuditTest` and
`GeeScaleAuditTest`; the positive-tail root/density/information is also frozen
by the R generator above.

1. Original-unit SVD at duplicated predictor scale 1e16 dropped an estimable
   intercept. Merely forcing scaled rank into that factor gave intercept2.58199
   instead of3; the repaired solver detects inaccurate projection/leverage and
   rejects explicitly. Narrow and wide rejection tests and valid original-norm
   controls are included. No silent wrong fit or alternate norm is returned.
2. CV with weights [1e200,1e-200,...] no longer rejects valid raw training
   weights after global-normalization underflow. Seed42 three-fold held-out risk
   is9 in the analytic dominant-weight limit. Numerical zero-mass folds do not
   contaminate aggregation; all training fits still must converge.
3. Poisson [1,2,3,4], intercept and all weights1e-20 now has SE3.16227766e9,
   not the floored-information SE1.58113883e7. The same true-information policy
   covers tiny binomial probabilities; no covariance uses stabilization floors.
4. All-zero Gaussian responses now have ML likelihood supremum+Infinity,
   matching OLS, rather than finite log density at a substituted variance.
5. An all-success group of1e22 trials with offset50 plus three ordinary
   failures now converges to beta.17001035679439 and LL-3.9724844032286.
   Predictor-aware p*q, stable response residuals, log1p deviance and
   complementary binomial density avoid subtracting an unrepresentable tail
   from one. Shared-family callers in Laplace/PQL were coordinated with their
   owner, who added separate integration regressions.
6. Gaussian GEE responses [1..6]*1e-12 now return Pearson phi3.5e-24 and naive
   covariance5.83333333e-25, not1e-12/1.66666667e-13. Scale-normalized score
   convergence handles the associated large raw score units. Exact zero
   residual variance uses a harmless common working scale and reports zero
   dispersion/covariance, rather than overflowing the sensitivity matrix.

The independent worker confirmed all six original gates pass after repair.
An intercept fit at binomial proportion1e-100 may still exhaust the default
100 iterations; it reports nonconvergence, and is never a valid timed success.

## Changed paths (owned scope)

- `src/main/java/org/jlinalg/internal/LeastSquaresSolver.java` (parent-authorized)
- `src/main/java/org/jlinalg/ols/Ols.java`, `OlsResult.java`
- `src/main/java/org/jlinalg/glm/Glm.java`, `GlmFamilies.java`, `GlmFamily.java`, `GlmResult.java`
- `src/main/java/org/jlinalg/gee/Gee.java`, `GeeConvergenceDiagnostics.java`
- `src/main/java/org/jlinalg/penalized/PenalizedRegression.java`, `PenalizedRegressionCrossValidation.java`
- `src/test/java/org/jlinalg/ols/OlsEdgeAuditTest.java`
- `src/test/java/org/jlinalg/glm/GlmEdgeAuditTest.java`, `FittingRAuditTest.java`
- `src/test/java/org/jlinalg/gee/GeeScaleAuditTest.java`
- `src/test/java/org/jlinalg/penalized/PenalizedFittingEdgeAuditTest.java`
- `src/test/resources/fitting-audit-v030/reference.properties`, `folds.txt`
- `src/benchmark/java/org/jlinalg/benchmark/FittingAuditBenchmark.java`
- `src/benchmark/r/fitting_audit_v030.R`
- `src/benchmark/resources/fitting-audit-v030/evidence.md`
- `docs/vignettes/linear-models-and-glms.md`

Independent read-only parent-core review and recheck are in ignored
`build/core-cross-review/REVIEW.md`; the parent fixed all five reported findings
and the original probes independently pass. No parent-owned core source was
edited by this worker. No commits, pushes or global Gradle builds were made.
