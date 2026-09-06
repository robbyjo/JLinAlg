# ZIP/ZINB outer-point caching benchmark

Measured 2026-09-05 against `8d25365f2ca405445fd0b01884155abd966f2aae`
(the previous adaptive-gradient implementation), on Windows, Intel Core
i9-14900K, Java 25, and JDistlib 0.10.1. No optimizer, evaluation budget,
convergence tolerance, likelihood approximation, or inference method changed.

## Results

The 644-observation Salamanders examples use a site random intercept (23
random coefficients). ZIP fits `count ~ mined + (1 | site)` with
`ziformula = ~ mined`; NB2 ZINB fits `count ~ spp + mined + (1 | site)` with
`ziformula = ~ spp + mined` and intercept-only dispersion.

| Model/path | Before median | After median | Speedup | Time reduction |
| --- | ---: | ---: | ---: | ---: |
| ZIP, full fit | 32.290 ms | 22.861 ms | 1.41x | 29.2% |
| ZINB, full fit | 211.470 ms | 115.625 ms | 1.83x | 45.3% |
| ZIP, prepared fit | 31.551 ms | 22.134 ms | 1.43x | 29.8% |
| ZINB, prepared fit | 211.161 ms | 115.305 ms | 1.83x | 45.4% |

Full fit includes symbolic preparation; prepared fit reuses it. Each median
pools 21 measured fits: three JVMs per revision, two warmups then seven
measurements per model/path per JVM. Revision order alternated
before/after, after/before, before/after; processes ran sequentially, without
profiling or competing test runs. Both versions used `-Xms512m -Xmx4g`, CPU
backend, AUTO optimization, and a maximum of four gradient workers (ZIP
remains serial under the existing small-parameter gate; ZINB uses four).
Inference was computed and printed separately, outside timed fits.

[Raw per-fit timings](benchmarks/zero-inflated-cache-2026-09-05.csv) are
checked in. These are machine/workload-specific macrobenchmarks, not a
universal speed guarantee or a statistical confidence interval.

The correlated-pedigree ZIP check (500 members, 1,500 rows, 1,000 random
coefficients) showed **no speedup**: medians of three fresh-process fits were
2.126370 s before and 2.140810 s after (0.7% slower). Both retained 4,900
equation/factor nonzeros, 327 objective evaluations, and exactly the same
reported likelihood, -2177.7034306278597. Sampled peak-heap deltas remained
approximately 310-312 MB. See the
[raw pedigree measurements](benchmarks/zero-inflated-pedigree-cache-2026-09-05.csv).
This change targets likelihood work, not sparse factorization.

## What changed and accuracy checks

An initial Flight Recorder profile attributed about 17% of sampled CPU time
to fixed-design dot products and 9% to log-gamma. Those calculations were
repeated during random-mode Newton iterations and line searches, even though
their inputs were constant at a given outer parameter point.

Each objective/gradient worker now owns reusable caches for count/zero fixed
predictors, NB2 sizes and their log/gamma terms, and structural-zero
probabilities/log probabilities for rows without zero-process random effects.
Caches refresh on **every** outer point, including finite-difference and
inference evaluations. Rows with zero-process random effects still recompute
their probabilities at each trial mode. The caches add eight row-length
double arrays per NB2 objective (five for ZIP), not dense random-effect
covariance storage. Existing arithmetic ordering is preserved.

All printed Salamanders coefficients, sizes, random-effect SDs, likelihoods,
standard errors, and evaluation counts were identical across the six JVM
runs. ZIP used 109 evaluations; ZINB used 1,131. Relative to a fresh glmmTMB
1.1.14 / R 4.6.1 run, absolute likelihood differences were 1.87e-8 (ZIP)
and 5.85e-8 (ZINB). This is agreement at the documented fitted solution,
not proof of a global optimum. Numerical standard errors are unchanged;
the pre-existing finite-difference inference method is not made exact by
this performance change.

The refreshed single-threaded R medians, excluding SE computation, were
30 ms for ZIP and 170 ms for ZINB (two warmups, seven measurements). R's
elapsed timer on this host was coarse; do not read sub-millisecond precision
into those numbers. The Java ZINB comparison uses four workers versus one
for R, so it is not an equal-core efficiency comparison. R medians including
SE computation were 40 ms and 250 ms; these are not compared against the
Java fit-only timings above.

`gradlew check` passed: 286 tests, zero failures/errors, three optional native
CHOLMOD tests skipped. Added coverage reconstructs an uncached block-diagonal
Laplace likelihood independently of the new caches, checks per-row dispersion
and offset predictions, compares serial/parallel gradients, and reuses
prepared structures with changing responses/designs, both with and without
zero-process random effects. Existing R/TMB and derivative checks also pass.

## Reproduce

From the repository root, with R/glmmTMB already installed:

```powershell
Rscript src/benchmark/r/salamanders_zero_inflated.R
.\gradlew.bat benchmarkSalamandersZeroInflated
.\gradlew.bat benchmarkZeroInflatedMixed
.\gradlew.bat check
```

The R script exports identical response/group/design matrices under
`build/tmp/salamander-comparison`; the Java task reads these, fails on
non-converged fits, and prints estimates plus individual timing samples.
Override `jlinalg.benchmark.salamanders.data`, `.threads`, `.warmups`, or
`.runs` with JVM system properties as needed.

For the before/after comparison, compile the baseline revision in a separate
checkout and run the **same current benchmark class** against either old or
new main classes. The measurement invocation was:

```powershell
java '-Djava.io.tmpdir=build/tmp' -Xms512m -Xmx4g -cp 'build/classes/java/benchmark;build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar' org.jlinalg.benchmark.SalamandersZeroInflatedBenchmark
```

Replace only the main-class directory in the classpath with the compiled
baseline classes. Do not time Gradle/JVM startup, R data export, or compilation
as part of the warmed Salamanders fits.
