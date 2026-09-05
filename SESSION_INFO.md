# Session information

Session date: 2026-09-05
Time zone: America/New_York (UTC-04:00)
Repository baseline: `main` at `9389698` (`origin/main`)

## Timestamped activity

- **2026-09-05 14:29 EDT / 18:29 UTC** — cloned the GitHub repository into
  the empty `E:\Projects\JLinAlg` workspace and reviewed the roadmap, survival,
  pedigree, sparse mixed-model, backend, test, and benchmark code.
- **2026-09-05 14:31 EDT / 18:31 UTC** — ran the pre-change targeted sparse Cox
  and sparse pedigree tests; the baseline passed.
- **2026-09-05 14:38 EDT / 18:38 UTC** — completed the first implementation and
  validation pass. Focused sparse Cox, dense Cox mixed, and sparse pedigree
  tests passed, including new repeated-incidence and dense/sparse
  conditional-mode parity cases. The updated benchmark source also compiled.
- **2026-09-05 14:39 EDT / 18:39 UTC** — synchronized the README, development
  inventory, Cox vignette, performance notes, and generated website page with
  the implemented behavior and remaining numerical boundary.
- **2026-09-05 14:43 EDT / 18:43 UTC** — completed final verification:
  `gradlew check` passed the full test and website suite, `benchmarkClasses`
  passed, and `git diff --check` found no whitespace errors. Three optional
  native-CHOLMOD tests were skipped because that native library is not packaged
  in this checkout.
- **2026-09-05 16:28 EDT / 20:28 UTC** — generated a deterministic
  800-row/200-member repeated-measures pedigree fixture and installed the R
  comparison packages into the ignored benchmark library.
- **2026-09-05 16:38 EDT / 20:38 UTC** — isolated R/JLinAlg parity with fixed
  Cox, shared-frailty, identity-pedigree, dense-Java, and sparse-Java controls.
  R's inverse relationship matrix matched the Henderson precision exactly.
- **2026-09-05 16:41 EDT / 20:41 UTC** — found that `coxmeMlist`'s sparse
  `Matrix` input path changed the structured-covariance fit. Supplying the same
  relationship as a base-R dense matrix reduced the maximum gene-beta
  difference to `1.76e-8`.
- **2026-09-05 16:46 EDT / 20:46 UTC** — completed the reproducible
  correctness/performance gate. All eight Java fits converged; maximum standard
  error difference was `1.02e-4`, and the warmed median was 0.035909 seconds
  for eight-worker JLinAlg versus 0.360 seconds for one-worker R (`10.03x`).
- **2026-09-05 16:50 EDT / 20:50 UTC** — reran the full `gradlew check`
  suite, benchmark compilation, deterministic fixture generation, both runtime
  benchmarks, the strict cross-runtime gate, and `git diff --check`; all passed.
  Three native-CHOLMOD-only tests remained skipped because the library is not
  packaged in this checkout.
- **2026-09-05 16:57 EDT / 20:57 UTC** — marked the sparse Cox mixed and
  pedigree TODO item completed, with approximation boundaries retained as
  explicitly separate follow-on enhancements.
- **2026-09-05 17:05 EDT / 21:05 UTC** — reviewed every remaining Major and
  Medium inventory item and separated stale documentation from substantive
  implementation gaps.
- **2026-09-05 17:14 EDT / 21:14 UTC** — added a reusable counting-process
  risk-set sweep, prepared Cox score scanning, cluster/relatedness variance,
  diagnostic residuals, proportional-hazards tests, and shared gamma frailty.
- **2026-09-05 17:18 EDT / 21:18 UTC** — replaced the sparse pedigree PEV
  cutoff with fixed-size inverse-column batches and exposed PEV/reliability plus
  selected-member lookup.
- **2026-09-05 17:21 EDT / 21:21 UTC** — completed hit-only xWAS robust
  follow-up, including a CLI output for RAPS, contamination-mixture, and
  PRESSO-style results.
- **2026-09-05 17:24 EDT / 21:24 UTC** — added the dedicated colocalization
  Markdown/website workflow and repaired stale MR, pedigree, Cox, and omics
  documentation.
- **2026-09-05 17:26 EDT / 21:26 UTC** — ran the Cox pipeline benchmark:
  28,748 prepared score tests/second versus 2,202 full refits/second, with a
  2,000-row start-stop fit taking 0.004636 seconds.
- **2026-09-05 17:27 EDT / 21:27 UTC** — ran the sparse pedigree uncertainty
  benchmark: all 500 PEVs were finite for 1,000 observations in 0.076384
  seconds using 32-column solve batches.
- **2026-09-05 17:28 EDT / 21:28 UTC** — verified fixed-theta gamma frailty
  coefficients, SEs, and modes and recurrent-event cluster sandwich SEs against
  R `survival` 3.8-3 under R 4.6.1; the checked-in tolerances are `1e-6`.
- **2026-09-05 17:34 EDT / 21:34 UTC** — ran the checked-in conventional R
  Cox scan comparator at 2,000 rows and 512 predictors. R processed 1,829
  predictors/second versus JLinAlg prepared score scanning at 28,748/second,
  a 15.72x measured advantage on this host.

## What changed

- Generalized `SparseCoxMixedModel` from a pedigree-only entry point to one
  unit-incidence random-effect term with any caller-supplied
  `SparsePrecisionMatrix`; the existing pedigree overload remains compatible.
- Added exact aggregation of the penalized score, dense cross-information, and
  random-information diagonal when multiple observations share one sparse
  coefficient. Unobserved pedigree ancestors remain in the sparse system.
- Added `CoxPedigreeFrailty.fitSparse(...)` as an explicit sparse facade. Its
  `PedigreeRandomEffectTerm` overload permits end-to-end sparse construction;
  the existing `fit(...)` remains the dense reference path.
- Added `CoxMixedSolver` plus result diagnostics for solver type, sparse
  coefficient count, sparse equation nonzeros, numeric-factor nonzeros, and
  lower-variance-bound singularity through `isSingular(tolerance)`.
- Extended `TopmedCoxBenchmark` to sample JVM peak heap and record sparse
  equation/factor diagnostics in its CSV artifacts.
- Added regression tests for repeated shared frailty, dense/sparse conditional
  mode agreement at effectively fixed variance, sparse diagnostics, and the
  repeated-individual pedigree facade.
- Added deterministic R fixture generation and an executable parity/speed gate
  against `coxme` 2.2-22. The R harness can select dense or sparse pedigree
  covariance representation and fixed variance components.
- Added a prepared fixed-variance dense reference path, warm-started prepared
  sparse scans (including fixed-variance scans), and routed the preferred
  sparse-Cox fallback to the portable CPU sparse kernel when native CHOLMOD is
  unavailable.
- Updated the benchmark to perform a complete untimed warm scan, report the
  actually selected per-fit backend, and compare the same fixed variance model
  in Java and R.
- Added `CoxCountingProcessPlan`, `FastCoxAssociation`, `CoxDiagnostics`, and
  `CoxGammaFrailty`, with streaming omics integration and reproducible
  throughput benchmarks.
- Added batched sparse diagonal PEV/reliability extraction with dense-reference
  and beyond-former-cutoff tests plus a scale benchmark.
- Added retained-hit-only xWAS instrument hand-off and optional robust CLI
  output, plus a standalone colocalization vignette in Markdown and HTML.

## Numerical contract and remaining work

The sparse Cox kernel does not materialize a dense random-effect covariance or
information matrix. It retains the exact penalized partial-likelihood score,
but uses the diagonal of the profiled random-effect information plus the full
sparse precision for Newton/Laplace calculations. Conditional modes therefore
target the dense score solution when both paths converge, while profiled
variance estimates and fixed-effect covariance are not promised to be
identical.

The sparse Cox mixed path still requires one-stratum right-censored data,
distinct event times, and one unit-valued incidence per observation. Exact
sparse Laplace information/determinants, multiple sparse terms, ties, strata,
start-stop sparse frailty, and a sparse-GRM convenience facade remain optional
follow-on enhancements outside the completed Major/Medium inventory. Fixed
Cox start-stop data use the new exact indexed sweep.

The historical TOPMed data were not present in this fresh checkout, so its
published large-cohort numbers were not regenerated. The new synthetic gate
was executed locally and is reproducible without those private inputs. Its
speed result is specific to the stated fixture, R's one-threaded scan, and
JLinAlg's eight independent scan workers.
