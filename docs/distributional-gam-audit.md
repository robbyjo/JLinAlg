# Distributional and GAM numerical audit (2026-09-08)

Scope: `gam` and `distributional`, including fixed/mixed beta and sparse
zero-inflated count models. Shared Laplace-engine changes were coordinated
with the GLMM owner. No new pedigree covariance optimizer or native-runtime
certification is claimed. The owned isolated suite passes 70/70 tests.

Post-release follow-up (2026-09-10) adds deterministic five-start outer-mode
search plus stationary-basin confirmation. The previously failing competing-
mode ZINB regression now converges and certifies the selected objective as best
and reproducible among those deterministic starts. This is a bounded
multi-start certificate, not a proof of the continuous global optimum or a
multi-mode random-effect integral. Observed-Hessian fixed-effect covariance is
reported conditional on nuisance variances at active bounds.

## Reproduced defects and repairs

| Case | Before | Repair / independent gate |
| --- | --- | --- |
| Weighted Gaussian GAM with offsets | fitted + residual differed from original response by minus offset | Restore offsets in returned fitted values; heteroskedastic smooth-fit regression preserves adjusted-response REML covariance/residuals |
| Beta density, y=mu=.4, precision=1e14 | LL 15.78515625 versus 15.912715295573715 | Stable JDistlib beta density in generic, specialized and mixed kernels |
| NB2 density, y=mu=2, size=1e14 | LL -1.41086363055 versus Poisson limit -1.30685281944 | Stable NB density and zero-mass log1p expression; 48 R beta/NB cases |
| Hurdle count information | NaN at mean 1000; artificial floor at mean 1e-10 | Zero-truncated-Poisson moments, expm1 and small-mean series; finite differences and independent R maximum |
| Hurdle/ZIP positive log PMF, y=mu=1e15 | Hurdle LL -16 instead of -18.411470281974225 | Stable Poisson density in fixed hurdle/ZIP and both sparse ZIP likelihood evaluators; R dpois gates through 1e15 |
| Beta/generic distributional maximumStep=1e-12 | False convergence after one iteration | Require scaled score stationarity as well as small step/likelihood change; maxiter=10 regression remains nonconverged |
| Beta random-intercept Laplace fit | LL 435.2985387 versus R 435.1762842886 | Marginal fixed-effect optimization plus observed random-mode Hessian; corrected LL 435.176284289528 |
| NB2 coefficient covariance | Stabilized block scoring metric reported as covariance | Invert observed Hessian including mean/size cross derivatives; full covariance max R error 1.72e-8 on 400-row variable-mean/size fixture |
| Sparse ZI stopping/profile checks | Small step or evaluation count could certify a nonstationary fit; unconverged inner mode could produce a successful profile | Central differences, projected-gradient checks, cold final stationarity audit, tighter inner modes, truthful evaluation counts; failed profile false/NaN; no covariance for uncertified fit |

The beta/NB density grid has maximum absolute R error 2.33e-9 on a large
negative log density; tests also enforce scaled errors. Extreme-density
agreement does not certify estimability or Hessian accuracy at arbitrarily
large precision. The NB optimizer uses a positive stabilized scoring metric
for steps, distinct from its final observed-information covariance.

The beta observed predictor information is Fisher information minus
`(1-2*mu)*score`; it can be negative away from the conditional mode. An
analytic/finite-difference grid explicitly tests that negative-curvature
case. The previous .122 likelihood discrepancy was a bug, not an acceptable
Laplace approximation difference. See [beta mixed models](beta-mixed-models.md).

## Accuracy-qualified warm timings

Five measured fits, Java CPU `-Xmx128m`, after three Java warmups. R mixed
fits use two warmups, R 4.6.1/glmmTMB 1.1.14. Every timed result must converge;
R mixed fits additionally require a positive-definite reported Hessian.
Checksums sum likelihood plus fitted fixed coefficients and are retained in
[`observed-run.txt`](../src/benchmark/resources/distributional-audit/observed-run.txt).

| Workload | Java median | R median | Accuracy |
| --- | ---: | ---: | --- |
| Fixed beta, variable precision, n=1500 | 4.663 ms | 50 ms | coefficient max error 6.06e-10 |
| Beta mixed, n=600, 40 groups | 244.961 ms | 200 ms | LL error 9.34e-10 |
| ZIP, n=450, 18 count and zero groups | 82.455 ms | 140 ms | LL error 2.05e-12 |
| ZINB, n=450, 18 count groups | 45.439 ms | 190 ms | LL error 1.01e-9 |

These are whole-call timings on the same data/likelihood, not equal-work
optimizer microbenchmarks: Java ZIP requests numerical inference, Java ZINB
is fit-only, and R includes sdreport. No general speed ratio is inferred.
The old 30,000-row beta mixed speed claims are withdrawn because they timed
the incorrect objective. Corrected large-pedigree throughput remains unmeasured.
Sparse mixed engines do not allocate an observation-by-observation covariance;
factor fill, random dimension and outer differentiation still limit scaling.
Weighted/generalized GAM paths remain dense, and no new scaling claim is made.

## Explicit unresolved ZINB counterexample

`SparseZeroInflatedMixedModelTest.outerCachesRefreshForDispersionOffsetsAndParallelGradients`
is the executable reproduction. It uses Java Random seed 7413922, 12 groups
of 20, x=-1+2*(row%20)/19, count mean exp(.45+.5*x+.75*sin(1.1*group)),
structural-zero probability .22, and NB size 2 (two-success Bernoulli sampler).
Count design is [1,x,sin(.73*row)], zero design [1,x]. The second refit adds
one to responses at row%7==0 and changes size design from [1,x] to
[1,cos(.17*row)]. Offsets start at .12*cos(.31*row) and are incremented by
.02 on each second pass; the test intentionally reuses its offset array
across the count-only and count-plus-zero-random-effect loops.

With count and zero grouped random intercepts, the changed-response refit
has competing conditional modes under small outer perturbations. Controls:
`ZeroInflatedMixedOptions(2000,100,1e-6,.4,1e-6,100,1e-4,1e4,20,null,BOUNDED_BFGS,threads)`
for threads 1 and 4. Both return false convergence, LL -394.345176925742,
880 objective calls. The regression checks identical serial/parallel fitted
state and call accounting and reconstructs the joint likelihood and observed
Hessian independently without caches. The other three cases still require
successful convergence; their LLs are -374.837034021561, -397.665820658931,
and -374.837046048419.

Mixture likelihood nonconcavity is a model property; inability to robustly
select/integrate competing modes is an unresolved implementation limitation.
This is not evidence of a certified boundary MLE. No MLE/inference or speed
claim is made for this failure. A separate profile regression restricts
inner mode iterations to one and verifies false convergence and NaN likelihood.

## Reproduction

From the repository root (R requires access to the existing user-library
dependencies on this host; no installation is needed):

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-distributional-audit-reference.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-distributional-poisson-extremes.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/distributional_audit_benchmark.R
.\gradlew.bat benchmarkDistributionalAudit
```

The Gradle task builds current sources with the pinned dependency. Full checks:

```powershell
.\gradlew.bat test --tests 'org.jlinalg.distributional.*' --tests 'org.jlinalg.gam.*'
.\gradlew.bat check benchmarkClasses assemble
```

The test suite uses checked-in fixtures; R is needed only to regenerate them.
Raw mixed R medians/likelihoods/checksums are also in
`src/benchmark/resources/distributional-audit/r-timings.properties`.
