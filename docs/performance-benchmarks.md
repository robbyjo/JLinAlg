# Performance benchmarks

`benchmarkAssociation` runs deterministic end-to-end macrobenchmarks for:

- covariate-reused block OLS;
- prepared-null block GLM score association;
- exact parallel per-predictor OLS refits.

`benchmarkMixedModels` measures sparse-equation REML and reports both equation
and factor nonzero counts. For workloads with at most 2,000 observations it
also runs the dense covariance reference path.

The task performs warm-up iterations, reports the median of measured runs, and
consumes a result coefficient to prevent dead-code elimination. It is intended
for workload and backend comparisons rather than nanosecond-scale JVM method
measurement.

```powershell
.\gradlew.bat benchmarkAssociation
.\gradlew.bat benchmarkMixedModels
```

Dimensions and execution can be controlled without editing source:

```powershell
.\gradlew.bat `
  "-Djlinalg.benchmark.rows=50000" `
  "-Djlinalg.benchmark.variables=10000" `
  "-Djlinalg.benchmark.covariates=10" `
  "-Djlinalg.benchmark.parallelism=16" `
  "-Djlinalg.benchmark.warmups=3" `
  "-Djlinalg.benchmark.measurements=7" `
  benchmarkAssociation
```

Output is CSV with median wall time and variables per second. Record the Java
version, CPU/GPU model, backend, thread environment, heap settings, and power
mode alongside published results. Run on an otherwise idle machine and compare
multiple repetitions before changing routing thresholds.

The default workload deliberately fits only 100 exact OLS models because that
path is a correctness/per-fit baseline; the fast paths scan the full requested
variable count.
