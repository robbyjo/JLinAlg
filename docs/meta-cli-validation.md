# Meta CLI validation and allocation profile

Run date: 2026-09-12. Java 25 on the Windows development host, portable CPU;
these are synthetic workload measurements, not end-to-end omics throughput.

The integrated `check javadoc executableJar` gate passed: 680 tests, 677
passed, zero failures, and three optional native CHOLMOD skips. Both commands
were also run from the built executable on the six-cohort tutorial files;
the output retained `++--?-` and the documented singleton/regression statuses.

Independent fixtures in `src/test/resources/meta/cli-reference.properties`
were generated with R 4.6.1 and metafor 5.0-1. CLI tests compare coefficients,
SEs, p-values, tau², and residual Q for fixed effects and REML, DL, and PM,
both with and without a moderator. PM's R root tolerance is explicitly set
to `1e-12`; its default tolerance was too coarse for this numerical comparison.
Other tests cover analytic equal-variance heterogeneity, missing cohorts,
singleton pass-through, per-row inference degrees of freedom, parallel parity,
unit changes, unsorted and compressed cohort files, a merge with more than
32 runs, duplicates, malformed inputs, protected outputs, and regression rank.

Reproduce the fixtures and checks from the repository root:

```powershell
Rscript src/test/resources/meta/generate-cli-reference.R
.\gradlew.bat test --tests org.jlinalg.cli.MetaCliTest --tests org.jlinalg.meta.MetaArrayTest
```

The committed `MetaArrayBenchmark` compares a per-feature list/scalar fit
with the missing-aware array batch. Inputs: 5,000 features, eight cohorts,
one cohort missing on half the rows. Both paths are single-threaded; input
generation is outside timing, while list construction or array preparation
is included. Each path is warmed once and measured three times. The table
reports medians, with results consumed and estimates compared.

| Model | Study-list time | Array-batch time | List allocated bytes | Batch allocated bytes |
| --- | ---: | ---: | ---: | ---: |
| Fixed | 7.912 ms | 1.960 ms | 24,940,272 | 2,620,480 |
| Random REML | 76.217 ms | 15.702 ms | 320,780,016 | 2,460,480 |

Maximum absolute pooled-beta differences were `1.39e-16` for fixed effects
and `7.32e-10` for REML. Allocation is measured with the JVM thread allocation
counter, not peak resident memory. The improvement combines avoiding study
objects with the specialized intercept-only arithmetic; it does not isolate
the cost of `List` itself. Meta-regression has an array interface but still
performs a general weighted matrix fit per feature, so these speed ratios
do not apply to regression. Sorting, disk I/O, JVM startup, and joining are
also outside this benchmark. Timings can vary on a shared host.

```powershell
.\gradlew.bat benchmarkClasses
java -Xmx512m -cp "build/classes/java/benchmark;build/classes/java/main;build/dependencies/jdistlib-all-0.10.2.jar" org.jlinalg.benchmark.MetaArrayBenchmark 5000
```
