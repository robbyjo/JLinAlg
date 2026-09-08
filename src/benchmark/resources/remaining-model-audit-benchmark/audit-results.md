# Survival, mediation and nonlinear audit (2026-09-08)

Scope: remaining code after dd45e80 in `org.jlinalg.survival`,
`org.jlinalg.mediation`, and `org.jlinalg.nonlinear`. No public method signatures
were changed. No global build configuration, site files, commits or pushes
were made by this audit worker.

## Reproduced defects and repairs

| Counterexample | Before | Repair / independent gate |
| --- | --- | --- |
| Nonlinear mean `f(b)=1e12*b`, y=[1,2,3], initial b=0 | b=0, SSE=14, converged=true | Scale-invariant projected Gauss–Newton check; b=2e-12, analytic SSE=2 |
| Nonlinear mean `f(b)=b*b`, y=[3,4,5], initial b=1, one iteration | Returned b=2.5 but SE=.846254 from the old Jacobian; iterations=2 | Recompute final Jacobian; analytic SE=.338502, iterations=1, converged=false |
| Mixed nonlinear `initialStep=.1`, one iteration, linear mean | Step controls were ignored and an undamped coefficient update used | Conditional-mean backtracking honors the controls; returns one-tenth of the complete update and converged=false |
| Mixed nonlinear assembly | Every fitted-value evaluation expanded sparse Z to an n-by-q dense array | Recover additive Zb from the linearized fitted vector; final linearization refreshed at returned parameters |
| Fixed Cox predictor multiplied by 1e-12 | Rescaled beta=0, converged=true; unscaled beta=.129902 | Internally center/scale predictors; restore beta, covariance and baseline parameterization |
| Disjoint start-stop risk sets, offsets [40,40,0,0] | Nonpositive denominator exception after large risks exited | Dynamically scaled active moments with cancellation-triggered rebuilding; log likelihood=-2 log(2), score=-1, information=.5 |
| A common Cox offset of 1e100 | Returned log likelihood 0 instead of -2 log(2) | Cancel the risk-set maximum in each log-likelihood contribution before summation |
| Sparse frailty: offset 40 on two rows censored before the first event | Nonfinite sparse-matrix exception | Reverse-time risk addition; matches the zero-offset conditional solution; an eta-range>300 locally scaled fallback also passes offset=1000 |
| Sobel a=1e160, b=1e-160, var(a)=1e300, var(b)=1e-300 | Infinite SE although the analytic SE is 1e10 | `hypot(b*sqrt(var(a)),a*sqrt(var(b)))`; invalid variances and unrepresentable outputs rejected |

Mixed stopping no longer interprets a worsening/rejected raw-SSE update as
convergence. Backtracking compares means holding the current BLUP fixed,
rather than comparing raw conditional SSE across different fitted covariance
models. A returned convergence flag also requires the variance fit to converge.

## Independent accuracy gates

Frozen inputs and references live in `src/test/resources/remaining-model-audit`.
The generator is `src/benchmark/r/survival_mediation_nonlinear_audit.R` and sets
the project R library before imports. Reference environment: R 4.6.1,
survival 3.8.6, lme4 2.0.6.

- Tied, stratified start-stop Cox: `survival::coxph` Efron and Breslow beta/SE/LL;
  largest beta difference 4.45e-16. Independent finite differences check score
  and information, and splitting each interval preserves the fit.
- Gamma frailty at theta=.4: independent R BFGS optimization of the penalized
  partial likelihood, including the gamma density and log-frailty Jacobian.
  Conditional mode error <=2.39e-7; partial likelihood and full Laplace
  normalization pass 1e-6 and 2e-6 tolerances, respectively. This compares the
  same Laplace approximation, not a differently integrated frailty estimator.
- Mediation: base-R `lm` and `lme4::lmer` for all three component models,
  comparing a, b, direct, total, indirect and Sobel SE. Maximum errors are
  1e-15 for OLS and 2.29e-9 for mixed models. Explicit 90% path intervals are
  tested with the reported finite denominator DF.
- Noisy exponential: `stats::nls` coefficients, final-Jacobian SEs and SSE.
  Maximum beta error 1.76e-10.
- Nonlinear pedigree ML: independent R dense `V=tau*Z*A*Z'+sigma*I`
  likelihood jointly optimized in nonlinear mean and log variances. The
  fixture contains founders, siblings and correlated parent-offspring effects.
  Maximum beta/log-variance error 2.26e-8; LL difference about 4e-14.

Fifteen new tests cover these references, adversarial units/offsets, analytic
counterexamples, finite differences, invalid inputs, and convergence budgets.
They are in `CoxAuditTest`, `MediationAuditTest`, and `NonlinearAuditTest`.
The final isolated JUnit package run against a frozen copy of the first
integrated main classes passed **35/35 tests, zero skips/failures** (12 test
containers). This includes the existing package regressions. The parent also
reported the first integrated 571-test run passing 568 with 3 optional skips.

## Warm wall-clock measurements

Ten warmups precede each measurement, followed by 50 fits (20 for the three
mixed/frailty cases). I/O, JVM/R startup and pedigree/design creation are outside
timing. Every Java mixed/nonlinear/Cox fit must converge; R optimization status
is checked, and NLS errors on nonconvergence. Each fit contributes coefficients
or effect estimates to a checksum. Identical frozen inputs are used. Java CPU
was run with `-Xms256m -Xmx1g`. Raw outputs are `measured-java.csv` and
`measured-r.csv` in this directory.
`measured-java-integrated.csv` records a second run against a frozen snapshot
of the parent's first integrated build: every accuracy error/checksum is
unchanged, while wall-clock variation illustrates the limits of short timings.

| Workload | n | Java ms/fit | R ms/fit |
| --- | ---: | ---: | ---: |
| Start-stop Cox, Breslow | 240 | .492 | 2.213 |
| Start-stop Cox, Efron | 240 | .346 | 2.139 |
| Gamma frailty, fixed variance | 36 | .253 | 15.535 |
| Ordinary mediation | 120 | .268 | .918 |
| Mixed mediation | 120 | 1.150 | 23.521 |
| Noisy nonlinear fixed fit | 240 | .320 | .977 |
| Nonlinear pedigree ML | 48 | 1.040 | 9.986 |

These are whole-fit API comparisons on this host, not universal kernel speedups.
R Cox/lm/lmer/nls include their formula/model-frame work; Java receives a
prepared numeric design/mean. The gamma reference is an independently written
R objective with numerical optimization/Hessian, and the pedigree reference
is a dense generic optimizer, not specialized `coxme`/`nlme` routines. In
particular, their ratios must not be advertised as specialized-package speed
advantages. Tests may be running concurrently elsewhere in the shared workspace.

## Reproduction

Run from `E:/Projects/JLinAlg` (R needs elevated access to its installed
user-library dependencies on this host):

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/survival_mediation_nonlinear_audit.R
$auditCp='build/remaining-audit-classes;build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar;'+((Get-ChildItem build/mixed-deps -Filter '*.jar').FullName -join ';')
New-Item -ItemType Directory -Force build/remaining-audit-classes | Out-Null
$auditSources=Get-ChildItem src/main/java/org/jlinalg/survival,src/main/java/org/jlinalg/mediation,src/main/java/org/jlinalg/nonlinear,src/test/java/org/jlinalg/survival,src/test/java/org/jlinalg/mediation,src/test/java/org/jlinalg/nonlinear -Filter '*.java'
javac --release 17 -Xlint:all -Werror -cp $auditCp -d build/remaining-audit-classes $auditSources.FullName src/benchmark/java/org/jlinalg/benchmark/RemainingModelsAuditBenchmark.java
java -cp $auditCp org.jlinalg.survival.CoxAuditTest
java -cp $auditCp org.jlinalg.mediation.MediationAuditTest
java -cp $auditCp org.jlinalg.nonlinear.NonlinearAuditTest
java -Xms256m -Xmx1g -cp $auditCp org.jlinalg.benchmark.RemainingModelsAuditBenchmark
```

`build/mixed-deps` is this workspace's cached JUnit classpath; on a clean checkout
use the dependencies resolved by Gradle. Full package tests can be reproduced
once shared sources are stable with `./gradlew test --tests 'org.jlinalg.survival.*'
--tests 'org.jlinalg.mediation.*' --tests 'org.jlinalg.nonlinear.*'`.

## Independent cross-review follow-up

Planck supplied a reproducible fixed-model counterexample with 100 rows,
`x=i%9-4`, `y=shift+2*x+.25*(-1)^i`, and mean `shift+b*x`, starting at `b=0`.
At `shift=1e15`, the former response-norm roundoff floor incorrectly certified
`b=0`, SSE `2714.25`, as converged after one iteration. The stationarity test now
has **no absolute response-scale floor**, and the fixed-model line search
requires strict SSE improvement. The repaired fit returns
`b=2.0014792899408285`, SSE `6.25`, after two iterations and honestly reports
nonconvergence because quantized mean evaluations still leave a nonzero
projected score. At `shift=1e12` it returns SSE `6.248516991735`, also with a
nonconverged precision plateau; at zero shift it converges normally. New
regressions cover those offsets and response units `1e-12`, `1`, and `1e12`
without enlarging tolerances. The existing noiseless exponential convergence
test still passes.

The second review counterexample used confidence level `Math.nextDown(1.0)`:
the Sobel CI upper endpoint incorrectly became infinite because the requested
CDF probability rounded to one. Both normal-product and mixed-path t intervals
now evaluate `(1-level)/2` directly in the upper tail. Independent R `qnorm`
and `qt` fixtures cover `.95`, `.999999999999`, and the largest double below
one. The example `a=1,b=2,varA=.04,varB=.09` now has finite upper endpoint
`6.146180537906798` at that last level. R fixture generation preserves 17 digits
so the level itself does not round to one on disk.

Additional reproduction command:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/nonlinear_mediation_crossreview.R
```

`measured-java-crossreview.csv` is a post-repair warm rerun of the benchmark
above (10 warmups; 50 or 20 fits). Accuracy and checksums are unchanged;
nonlinear NLS took `0.310804` ms/fit and pedigree nonlinear `1.288320` ms/fit.
As in the original measurements, wall-clock variation is not a kernel speed
claim and all fitted benchmark cases must actually converge.

Final isolated compilation passed `javac --release 17 -Xlint:all -Werror`.
The post-cross-review JUnit run passed **38/38 owned-package tests**, with
zero skips/failures (105.773 seconds; includes existing backend checks).
`test-summary.txt` records the classpath and result. Planck independently
reran both repaired counterexamples and confirmed them resolved, also
retaining passing Cox risk/nuisance, nonlinear mixed returned-state, and
Sobel reciprocal-unit probes.

## Limits retained

- Sparse Cox supports one stratum, right censoring, distinct event times and
  unit incidence. Its exact penalized score uses a diagonal random-information
  approximation plus sparse precision; Laplace variances/covariance are not
  claimed equal to dense full-information fits.
- Counting-process rebuilding can be quadratic for adversarial active-set
  churn. Sparse extreme-predictor fallback is O(n*events*p). Neither constructs
  an n-by-n matrix. A very large uncentered baseline hazard may still be outside
  floating-point range; partial-likelihood invariance does not guarantee a
  representable baseline on every original offset scale.
- Nonlinear random effects are additive on the response scale. Nonlinear
  parameter random effects, non-Gaussian integration, and general `nlme` parity
  are outside this path. REML uses the local first-order linearization; the
  independent joint nonlinear pedigree reference specifically uses ML.
- Sobel inference assumes zero a/b cross-model covariance and uses a normal
  approximation for the product. There is no bootstrap, robust cross-equation
  covariance, correlated cross-equation random-effect fit, or random-slope
  product mediation. Undefined zero-SE/zero-effect tests remain NaN.
