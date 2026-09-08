# Remaining GLMM, pedigree, and REML audit

2026-09-08, post-dd45e80 shared workspace. Java 25 targeting Java 17; R 4.6.1,
lme4 2.0.6. No full-project build, commit, or push was performed by this owner.
The previously validated adaptive quadrature implementation was not modified.

## Reproduced defects and repairs

1. **P1, marginal Laplace optimization.** Both old dense and sparse engines
   estimated beta inside joint penalized conditional-mode equations, then varied
   only variance/family parameters. On `rare.tsv` (384 rows, 32 groups, nonzero
   offsets), both reported convergence at beta `[-2.6930514912,.7809940616]`,
   variance `1.2319302917`, LL `-106.7587015664`. The independently optimized
   first-order Laplace objective instead has beta
   `[-2.9662415325,.8162423548]`, variance `1.3244160264`, LL
   `-106.4090535610`. Bounded BFGS now optimizes beta, log variances, and supported
   family parameters jointly; inner random modes hold beta fixed. Predictors
   are scaled, tiny variance starts get a finite-scale restart, and a small
   projected marginal score is required. No success is inferred just from a
   small likelihood change or exhausted line search.
2. **P1, observed curvature and inference.** The determinant previously used
   Fisher/stale working curvature and beta covariance came from the working
   equations, omitting full marginal/nuisance derivatives. The final accepted
   random mode now supplies observed curvature, with likelihood step-halving
   during mode search. `LaplaceFamilyDerivatives` permits exact score/observed
   information; other noncanonical families numerically differentiate their
   working score. Interior beta covariance comes from the full marginal Hessian.
   Invalid or nonconverged inference is entirely NaN. Dense covariance bases
   are factorized once and use the same corrected coefficient-space engine.
3. **P2, ancestry cycles in direct sparse construction.** Entries a<-b and
   b<-a were accepted by `ofUninbred`, giving the positive-definite matrix
   `[[5/3,-4/3],[-4/3,5/3]]`. Positive definiteness does not make that a pedigree.
   Both direct sparse paths now validate the directed ancestry graph in linear
   time without constructing dense A.
4. **P2, artificial REML identification.** Two identity bases with different
   labels both received variance 5/3 and `converged=true`; a ridge also enabled
   finite Satterthwaite inference. One all-observations random-intercept basis
   alongside an unrestricted intercept and residual covariance likewise reported
   an arbitrary genetic variance as converged. Validation now checks normalized
   covariance-basis rank on the REML contrast space, plus normalized information
   rank before inversion. ML and REML are checked on their distinct spaces.
   Duplicate labels and nonidentifiable variance directions are rejected.
5. **P2, indefinite covariance bases.** With y=(-1,1,-2,2), intercept-only X,
   known covariance 10I, and estimated basis diag(-1,1,1,1), the old engine
   returned `converged=true`, variance 1.23409804e-5. Each immutable basis and
   known contribution is now checked for positive semidefiniteness, independently
   of whether their sum happens to be positive definite. Validation is cached
   for reused bases; diagonal checks avoid unnecessary eigendecomposition.

Final bound-control probe also reproduced finite-difference evaluations outside
a valid nuisance interval [0,1e-6]. Gradient and information stencil widths now
respect both bounds; analytic regression tests cover a boundary score and a
narrow interior Hessian.

The subsequent predictor-aware GLM family integration exposed a caller defect:
at eta50 with an all-success weight1e22, old PQL produced z50 and covariance
1e-150 instead of z51 and .5184705528587072. Laplace's mean-only density failed
the finite-MLE positive-tail example (beta -13.60, LL -2.22e6, nonconverged).
PQL now uses family working-response/precision and predictor-aware deviance;
Laplace uses predictor-aware likelihood. Neither clips working precision to an
arbitrary floor; nonrepresentable working models are explicitly rejected.
With near-zero bounded random variance, the tail fit now converges at beta
.17001035669029554, LL -3.972484403228629, against independent R values
.17001035679438581 and -3.9724844032274467. New tail fixture keys and two
regressions cover this interaction without claiming a general exact GLMM fit.

**Independent cross-review P1, duplicate random covariance:** two differently
named grouped terms with the same grouping could return `converged=true` with
arbitrary variance splits and identical likelihood. The new
`CovarianceIdentification` preparation gate applies `Z Q^-1 Z'` through sparse
solves, reusing the precision factorization, and checks normalized action rank
by two-pass modified Gram-Schmidt. Four fixed-seed continuous probes are used
for large observation counts; n<=8 uses the complete coordinate basis. No n*n
matrix is created; auxiliary observation storage is O(n times number of terms).
Actions of linearly dependent covariance operators cannot have independent
rank. An inconclusive/near-dependent numerical result is conservatively rejected
rather than labeled an exactly proved algebraic dependence. A single nonzero
design with SPD precision needs no probing. This checks structural variance
bases, not general likelihood/nuisance identification for every family.

The original 16-count public probe now rejects duplicate terms for both starts
[1,1] and [.1,2], while its single-term fit remains converged with variance
1.6290117990 and LL -30.46986176789676. Four regression methods cover those
cases, equivalent column permutation/precision rescaling, a three-term
K3=K1+K2 dependency, independent controls, and 100,000-row sparse covariance
actions. Ampere independently recompiled and reran the original unchanged
probe: both duplicates reject, single-term/crossed controls pass, and permuted,
rescaled and three-term dependent covariance controls reject.

The final parent spot-check also reproduced an overflowing covariance scale
product: multiplying two 1e160 column scales erased a representable subnormal
slope variance. Exponent-separated covariance back-transformation now avoids
that product (and reciprocal overflow for tiny scales). A full rare-binomial
fit with slope units1e160 failed the new regression before the fix and passes
after it; direct tests additionally cover1e-160 and mixed1e-300/1e300 scales.

Ancillary corrections: scale-relative dense eigenvalue rank tolerance; explicit
rejection of families with unmodeled estimated residual dispersion; tunable
family fits sharing mutable family parameters are serialized. The old
`initialLogVarianceStep` control now initializes the nuisance-coordinate BFGS
metric rather than silently becoming unused.

## Independent accuracy evidence

`generate-remaining-mixed-reference.R` writes all fixtures under
`src/test/resources/r-reference/remaining-mixed/`:

* Rare-binomial and Poisson grouped fits: independently find every scalar random
  mode with `uniroot`, evaluate the first-order Laplace determinant, and optimize
  the complete objective in base R. A numerical joint Hessian supplies covariance.
  `glmer(nAGQ=1,tolPwrss=1e-12)` independently agrees: LL differences 1.38e-11
  and 6.14e-11. Default R mode tolerance gave a small discrepancy, so both sides
  were tightened instead of conflating it with approximation error.
* Unbalanced Gaussian ML and REML: dense R covariance/GLS/log-determinant
  optimization, separately cross-checked with `lmer`. Includes coefficients,
  variance components, full fixed covariance and full likelihood constants.
* Fixed-dispersion PQL: independently relinearized dense R REML with residual
  covariance diag(1/workingWeight), without a free residual scale. Coefficient,
  random predictor, and variance comparisons pass 2e-6. This is not claimed as
  parity with MASS::glmmPQL's additional residual-scale estimation or with ML.
* Inbred pedigree, one known parent, shuffled input order and repeated phenotypes:
  R tabular A and `solve(A)` match Java sparse precision within 2e-14. Dense and
  sparse Laplace paths agree with unobserved ancestors retained, including
  observation-level random predictors and variance.

20 new regression methods cover those references and the reproduced failures,
large design units (1e12), tiny variance starts (1e-10), mode/outer exhaustion,
unavailable inference, and invalid dispersion. **66/66 isolated tests pass**:
63 in GLMM/pedigree/REML plus three independently rerun beta mixed reference tests.
Production and tests compile with `--release 17 -Xlint:all -Werror`.

Coordination with the distributional owner supplied exact beta derivatives.
The original 600-row/40-group glmmTMB beta fixture now returns LL
435.176284289528, an error of 9.31e-10, versus the previously tolerated ~0.122.
That old discrepancy was a conditional-beta/Fisher-curvature implementation
defect, **not** first-order Laplace approximation bias. Beta/ZI source repairs,
their own reference fixtures, and their timing reports remain distributional-owned.

## Warmed performance, same accuracy gates

Three warmups, seven recorded runs, full fitting and fixed-effect inference.
Every Java and R fit must converge and pass beta/variance/covariance errors
2e-6 to 3e-6 and full LL error 2e-8; no timings from failed fits are accepted.
Java uses nanosecond wall time; short R Gaussian fits use batches of 20 and
report seconds per fit, avoiding misleading zero timings from R's clock.
All timed results contribute to checksums. Input/basis creation is outside the
Java timer; dense Laplace factorization and sparse numeric preparation are inside.
R includes formula/model setup. Reused REML covariance validation is cached.

| Case | Java sparse seconds | Java dense seconds | R lme4 seconds |
| --- | ---: | ---: | ---: |
| Rare Bernoulli Laplace | 0.016597 | 0.019843 | 0.060000 |
| Poisson Laplace | 0.033664 | 0.035692 | 0.060000 |
| Unbalanced Gaussian REML | — | 0.003006 | 0.007000 |
| Unbalanced Gaussian ML | — | 0.002820 | 0.007000 |

Java full LL errors in these refreshed runs are below 7e-13. Raw runs and checksums are
in `java-warm.tsv` and `r-warm.txt`. The R file preserves harmless host locale
startup warnings. Shared-host jitter is visible; these are small-workload
measurements, not universal speedups. No timing claim is made here for PQL,
large pedigrees, or distributional beta/ZI models. Prior timings from the wrong
marginal objective are not used as an accuracy-equivalent performance baseline.

## Reproduce from the repository root

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-remaining-mixed-reference.R
./gradlew benchmarkRemainingMixedAudit
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/remaining_mixed_benchmark.R
```

The R scripts prefix `.libPaths(c('build/r-library',.libPaths()))`; on this host
they require elevated read access to already installed user-library dependencies.
No package installation is needed. The Gradle task is registered by the parent
for the integrated build. Regression checks run with `./gradlew check` after
all owners stabilize. The saved audit used isolated compilation and JUnit
discovery; no global build was attempted during concurrent edits.

## Changed paths and remaining limits

Production: `glmm/{GlmmPql,GlmmLaplace,SparseGlmmLaplace,GlmmLaplaceOptions,GlmmLaplaceResult,
LaplaceOptimization,LaplaceFamilyDerivatives,CovarianceIdentification}.java`,
`pedigree/PedigreeRandomEffectTerm.java`, `reml/{Reml,VarianceComponent}.java`.
Tests: `RemainingLaplaceAccuracyTest`, `LaplaceIdentificationTest`, `RemainingPedigreeAccuracyTest`,
`RemainingRemlAccuracyTest` in their respective owned packages.
Evidence: the R generator/directory above, `RemainingMixedBenchmark.java`,
`src/benchmark/r/remaining_mixed_benchmark.R`, and this report/raw output directory.
Source vignette: `docs/vignettes/pedigree-and-glmm.md` (general Laplace/pedigree/REML
paragraphs only; preserve the distributional owner's ZI section).

PQL remains a working-likelihood approximation; first-order Laplace remains
approximate and can be biased for sparse binary groups. No general pedigree
adaptive quadrature, global-optimum certificate, exact zero-variance Laplace
fit, or finite-sample Wald coverage is claimed. Inference at configured bounds
conditions on those bounds. Generic family derivatives must match the supplied
likelihood. The sparse engine retains sparse precision/design, but dense
`Pedigree.of` still constructs A and dense REML remains observation-quadratic
in storage. Direct sparse pedigree construction trusts supplied inbreeding
coefficients (or the uninbred assertion); it validates ancestry but does not
independently derive biological inbreeding consistency.
