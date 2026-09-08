# JDistlib 0.10.2 dependency-upgrade validation

Date: September 8, 2026. Baseline: JLinAlg commit
`6cc2bb138c80d02018c199e161cfbe519861af4a` with JDistlib 0.10.1.
Current source builds pin [JDistlib 0.10.2](https://github.com/robbyjo/JDistlib/releases/tag/v0.10.2)
and verify `jdistlib-all-0.10.2.jar` against SHA-256
`8b50aaa11d2cbd512525e0d74f8e609f40db9f4fc179ccb7a885a82257e0ac81`.
The downloaded artifact matched. Published JLinAlg v0.3.0 binaries and its tag
are unchanged and still contain/use JDistlib 0.10.1; build current source to
use this upgrade.

## Outcome

No substantial accuracy regression was found. The original 609-test baseline
had 606 passes, three optional native CHOLMOD skips, and no failures. The
upgraded build has 610 tests: 607 passes, the same three skips, and no failures.
All nine accuracy-gated audit benchmarks and library/CLI/source/Javadoc
assembly passed. No production-algorithm repair was necessary.

One historical R-reference assertion initially failed because **JDistlib became
more accurate**, not less. We retained the original R data, corrected that
assertion to an independent finite-product reference, and added a 45-case
regression grid. No existing numerical tolerance was widened.

## Negative-binomial difference: R is not always the exact reference

For count `y=2`, mean `mu=100000`, and size `s=1e14`:

| Calculation | Log density |
|---|---:|
| JDistlib 0.10.1 / R 4.6.1 `dnbinom` | -99977.66729625262 |
| JDistlib 0.10.2 | -99977.66724625262 |
| Independent finite product | -99977.66724625262 |

For integer y the exact log density can be written as

```text
sum(k=0..y-1, log1p(k/s)) + y*log(mu) - log(y!)
    - (s+y)*log1p(mu/s).
```

The [R small-count implementation](https://raw.githubusercontent.com/wch/r-source/trunk/src/nmath/dnbinom.c)
uses a `-mu` approximation in its `y < 1e-10*s` branch. At the example above,
the omitted correction is about `0.00005` in log density. The new JDistlib
value agrees with the independent formula, evaluated separately in R.

The new grid crosses counts 0, 1, 2, 10, 100; means 0.01, 2, 100000; and sizes
1e8, 1e12, 1e14. Maximum scaled absolute error, defined as
`abs(actual-reference)/max(1,abs(reference))`, improved from
`5.0050061e-8` to `1.5007571e-9`. The largest absolute R-versus-finite-product
discrepancy on this grid is approximately `0.005`.

Residual limitation: outside the improved small-count branch, both dependency
versions retain up to `1.51e-9` scaled error on these cases. This is not a new
regression and does not warrant a local likelihood replacement for this quick
upgrade. The added test uses `5e-14` scaled tolerance in the improved branch
and an explicit `2e-9` elsewhere. This is a bounded fixture check, not a
claim of uniform tail accuracy over the parameter space.

The raw [grid](../src/test/resources/r-reference/jdistlib-upgrade-nb.tsv)
retains both R and finite-product values. Its
[generator](../src/test/resources/r-reference/generate-jdistlib-upgrade-reference.R)
does not call JDistlib to construct expected values. The original 48-case
beta/NB grid now has maximum absolute error `2.33e-9` against its R values
with the single independently justified correction above.

## Model-level accuracy

The full suite covers JDistlib-dependent distributions, inference, matrix
solves, OLS/GLM/GEE, mixed and pedigree models, beta and zero-inflated models,
SEM, meta-analysis, time series, MR/genetics, survival, mediation, nonlinear
and penalized/nonparametric regression. Existing independent R/reference
fixtures were held fixed except for the explicitly corrected assertion.
This is downstream integration coverage, not an exhaustive JDistlib audit.

Representative benchmark errors were unchanged at printed precision:

| Check | Error against existing independent reference, both versions |
|---|---:|
| Fixed beta coefficients | 6.06e-10 maximum absolute |
| Beta mixed log likelihood | 9.34e-10 absolute |
| ZIP / ZINB log likelihood | 2.05e-12 / 1.01e-9 absolute |
| Cox Breslow / Efron coefficients | 4.44e-16 / 1.67e-16 maximum absolute |
| Gaussian / mixed mediation | 9.99e-16 / 2.28e-9 maximum reference error |
| Nonlinear / nonlinear pedigree | 1.76e-10 / 2.26e-8 maximum reference error |
| Coloc / single-effect regression | 5.55e-17 / 8.88e-16 maximum scaled absolute |
| Winner's-curse correction | 2.17e-12 maximum scaled absolute |
| LD41 | 4.22e-15 relative |

Printed checksums also matched for OLS, Poisson, GEE, LASSO path/CV, the
distributional and genetic benchmarks, scan/BH pipelines, and the rank-two
and rank-three kernel tests. Matching checksums alone are not proof of every
output: the separate regression assertions remain the accuracy gates.

## Quick timing observation

Same Windows host, Java 25 (Java 17 release target), 32 logical processors,
sequential Gradle tasks, existing benchmark warmups/repetitions. These are
one before/after run, not interleaved controlled performance experiments.

| Workload | 0.10.1 median ms | 0.10.2 median ms |
|---|---:|---:|
| OLS, n=3000 | 0.368380 | 0.347720 |
| Poisson GLM, n=3000 | 1.603830 | 1.587285 |
| GEE fixed, n=3000 | 4.701820 | 4.679250 |
| LASSO CV, n=3000 | 0.855145 | 0.910300 |
| Fixed beta, n=1500 | 6.361 | 4.754 |
| Beta mixed, n=600 | 249.692 | 188.816 |
| ZIP, n=450 | 81.830 | 77.152 |
| ZINB, n=450 | 45.493 | 44.370 |
| Cox Breslow | 0.460338 | 0.547170 |
| Mixed mediation | 1.283815 | 1.439235 |

Several workloads improved, especially the beta fits in this observation;
others were slower. Do not infer universal speedups or attribute small changes
to the dependency without repeated measurements. R timings printed by the
distributional benchmark are **stored prior R measurements**, not freshly
timed R runs in this upgrade. The fresh R execution here generated the
independent NB grid. The [historical release report](release-0.3.0-audit.md)
retains the full R-comparison provenance and estimator limitations.

Raw Gradle benchmark logs: [before](../src/benchmark/resources/jdistlib-upgrade/before.txt)
and [after](../src/benchmark/resources/jdistlib-upgrade/after.txt).
The before log predates the extra finite-product test; the after log follows
the justified assertion correction. Java 25 emitted dependency native-access
and deprecated Unsafe warnings on stderr; no native/GPU certification is
implied by the portable success.

## Reproduce

From the repository root, with Java 17+ and R available:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-jdistlib-upgrade-reference.R
.\gradlew.bat check assemble benchmarkCoreStatisticsAudit benchmarkFittingAudit benchmarkDistributionalAudit benchmarkRemainingMixedAudit benchmarkRemainingModelsAudit benchmarkGeneticAudit benchmarkPipelineOlsAudit benchmarkPipelineBhAudit benchmarkPipelineKernelAudit --no-parallel --max-workers=1
```

The original baseline can be reproduced in a separate checkout of commit
`6cc2bb1`; do not reset a checkout containing your own work. Current reproduction
snippets use Gradle's pinned classpath rather than naming an obsolete JAR.
