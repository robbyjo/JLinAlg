# Formulas, model tables, and compute backends

## CLI-only formula and backend preflight

```bash
java -jar jlinalg-0.3.4.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "BMI ~ Sex * treatment + Age + <omics> + (1|site)" --backend preferred --threads 96 --block-size auto --out bmi-model.csv --dry-run

java -jar jlinalg-0.3.4.jar --pheno phenotype.csv --omics expression.csv --id SampleName --formula "BMI ~ Sex * treatment + Age + <omics> + (1|site)" --backend cpu --threads 96 --block-size auto --out bmi-model.csv --overwrite
```

`--dry-run` reports alignment, resolved model, exact-refit versus null-model
routing, block size, worker capacity, and backend without fitting. Strict
backend selections fail when unavailable rather than silently falling back.
The default `preferred` policy records the concrete backend selected.

## Compile a fixed-effect formula once

`ModelTable` is a lightweight column store used only while compiling a model.
The compiled object owns primitive contiguous arrays for fitting:

```java
ModelTable table = ModelTable.builder(6)
    .numeric("y", 11, 13, 15, 17, 19, 21)
    .numeric("x", 0, 1, 2, 3, 4, 5)
    .numeric("logExposure", 0, 0, 0, 0, 0, 0)
    .numeric("w", 1, 2, 1, 2, 1, 2)
    .categorical("sex", "F", "M", "F", "M", "F", "M")
    .build();

CompiledFormula model = Formula.compile(
    "y ~ x * sex + offset(logExposure)", table,
    new FormulaOptions(ContrastCoding.TREATMENT, "w"));

System.out.println(model.coefficientNames());
OlsResult fit = model.fitOls(
    OlsOptions.defaults(), BackendPolicy.PREFERRED);
```

`x * sex` expands main effects plus the interaction. Treatment contrasts use
the first encountered factor level as reference. Select `ContrastCoding.SUM`
for sum-to-zero contrasts. Compile outside repeated loops so formula parsing,
factor discovery, contrast expansion, missing-row compaction, and allocation
are not on the numerical hot path.

The same compiled fixed design can fit a GLM with `fitGlm` and an explicit
family/options/backend.

## Mixed formulas

Random-intercept, slope, independent double-bar, nesting, and correlated block
syntax follow familiar lme4 conventions:

```java
ModelTable repeated = ModelTable.builder(8)
    .numeric("y", 1, 2, 3, 4, 5, 6, 7, 8)
    .numeric("time", 0, 1, 0, 1, 0, 1, 0, 1)
    .categorical("site", "a", "a", "a", "a", "b", "b", "b", "b")
    .categorical("subject", "u", "u", "v", "v", "u", "u", "v", "v")
    .build();

CompiledMixedFormula independent = MixedFormula.compile(
    "y ~ time + (1 + time || site) + (1 | subject)", repeated);
SparseLinearMixedModelResult sparse = independent.fitSparse(
    RemlOptions.builder().initialVariances(1, 1, 1, 1).build(),
    BackendPolicy.PREFERRED);
```

Each `||` coefficient gets an independent variance term. A separate formula
such as `y ~ time + (1 | site/subject)` expands `site` and `site:subject`;
ensure the number of initial variances equals the actual compiled random terms
plus one residual term.

Use a single bar for an estimated within-group covariance:

```java
CompiledMixedFormula correlated = MixedFormula.compile(
    "y ~ time + (1 + time | subject)", repeated);
CorrelatedLinearMixedModelResult fit = correlated.fitCorrelated(
    RemlOptions.defaults(), BackendPolicy.PREFERRED);
```

`fitSparse` handles independent terms. `fitSparseUnstructured` and the
`fitCorrelated` convenience view use the sparse grouped-covariance optimizer
for correlated blocks.

Map an independent random-intercept group to pedigree precision at compile
time. Unobserved pedigree members remain random coefficients, and ordinary
independent terms may be fitted alongside the mapped term:

```java
CompiledMixedFormula pedigreeModel = MixedFormula.compile(
    "trait ~ age + (1 | individual) + (1 | site)", table,
    MixedFormulaOptions.builder()
        .pedigree("individual", pedigree)
        .missingDataPolicy(MissingDataPolicy.OMIT)
        .build());
SparseLinearMixedModelResult pedigreeFit = pedigreeModel.fitSparse(
    RemlOptions.defaults(), BackendPolicy.PREFERRED);
System.out.println(Arrays.toString(pedigreeModel.retainedRows()));
```

`MissingDataPolicy.OMIT` computes one mask across the response, every fixed
and random predictor, grouping columns, weights, and offset. This keeps all
compiled matrices and pedigree incidence rows aligned. The default
`MissingDataPolicy.ERROR` reports the first incomplete model row. Pedigree
mappings currently accept `(1|group)`. Mapped pedigree terms may be combined
with an unstructured grouped block; their covariance scales are optimized
jointly.

For an unstructured block with more than two coefficients,
`profileCorrelation(block, firstEffect, secondEffect, ...)` profiles any
selected effect pair. The implementation preserves the full covariance family
by permuting that pair into the leading Cholesky coordinates, including the
plus/minus-one correlation boundaries.

## Backend policy

`BackendPolicy.PREFERRED` follows this order:

1. GPU with JDistlib automatic workload routing;
2. oneMKL;
3. OpenBLAS;
4. portable Java CPU.

Pass it explicitly when reproducibility metadata matters:

```java
OlsResult fit = Ols.fit(
    y, x, OlsOptions.defaults(), BackendPolicy.PREFERRED);
System.out.println(fit.backend().selectedBackend());
System.out.println(fit.backend().deviceDescription());
```

For deterministic provider-specific checks use `BackendPolicy.CPU`. Strict
`GPU`, `CUDA`, `OPENCL`, `VULKAN`, `ONEMKL`, and `OPENBLAS` policies fail if the
requested provider is unavailable; they do not silently fall back.

Do not multiply outer Java threads by native BLAS/GPU worker pools without
benchmarking. Association defaults use one submitting thread with accelerator
routing, while `AssociationEngineOptions.cpuParallel()` deliberately selects
portable CPU for safe outer parallelism.

## Benchmark on the deployment machine

```powershell
.\gradlew.bat benchmarkAssociation
.\gradlew.bat benchmarkMixedModels
```

Tune benchmark sizes through the documented `jlinalg.benchmark.*` system
properties. Compare warm runs, record the concrete provider/device, and keep
the scientific model and numerical tolerances identical across backends. See
the [benchmark guide](../performance-benchmarks.md) for the complete protocol.
