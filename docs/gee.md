# Generalized estimating equations

JLinAlg's GEE implementation fits marginal generalized linear models to clustered or repeated observations. It uses contiguous row-major cluster blocks internally, JDistlib factorization backends, and cluster-level accumulation, so memory grows linearly with the number of observations rather than with a dense whole-study covariance matrix.

## Basic fit

```java
GeeOptions options = GeeOptions.builder()
    .correlation(GeeCorrelation.EXCHANGEABLE)
    .covariance(GeeCovariance.ROBUST)
    .build();

GeeResult fit = Gee.fit(
    response, design, subjectId, visit,
    GlmFamilies.binomial(), weights, offset,
    options, BackendPolicy.PREFERRED);
```

`design` includes the intercept when one is wanted. Cluster identifiers need not be sorted; `visit` supplies the within-cluster wave and must be unique inside a cluster. The overload accepting a flat row-major design avoids allocation of a `double[][]`. `GeeFormula.fit` accepts a reusable `CompiledFormula`, or declarations can be embedded directly as `y ~ treatment + visit + cluster(subject) + wave(visitNumber)`.

The result exposes coefficients, selected covariance and Wald inference, all covariance variants, fitted means, response/Pearson/deviance/working/standardized residuals, working-association and scale estimates, family-specific QIC/QICu/CIC/QICC, retained-row mapping, cluster influence and deletion diagnostics, convergence norms, new-data prediction, linear-contrast tests, tidy coefficient rows, JSON serialization through `toJson()`, and backend provenance.

## First-release surface

- GLM families and links through `GlmFamily`, with observation weights, offsets, estimated or fixed dispersion, missing-row omission, and CPU/accelerated backend selection.
- Independence, exchangeable, AR(1), m-dependent, Toeplitz, unstructured, fixed, and user-designed working association structures.
- Correlation and binary alternating-odds-ratio association models.
- Model-based and robust sandwich inference plus QIC-family model-selection criteria.
- Rectangular and contiguous row-major entry points, formula integration, immutable defensive-copy results, and reproducible R-reference tests.

## Extended surface

The five follow-on areas are included:

1. Small-sample inference: cluster degrees-of-freedom scaling, Mancl-DeRouen, Kauermann-Carroll, Fay-Graubard, and exact delete-one-cluster jackknife covariance; `CLUSTER_T` selects coefficient t and joint F inference with `clusters - parameters` denominator degrees of freedom.
2. Bias reduction: iterative bias-reduced estimating equations, a one-step bias correction, Jeffreys-penalized updates, and one-step/hybrid Jeffreys modes.
3. Nuisance models: observation-level dispersion regression with configurable identity/log/logit links, fixed association matrices, user-defined association design matrices, configurable identity/Fisher-z association links, constrained odds-ratio links, and damped alternating updates.
4. Categorical outcomes: `OrdinalGee.fit` supports proportional odds with working correlation or local odds ratios, `fitPartial` supports partial proportional odds, and `NominalGee.fit` provides marginal multinomial logits with arbitrary cluster dependence protected by a working-independence sandwich.
5. Clinical and diagnostic reporting: `ClinicalGee` provides reference-grid means, prediction, delta-method uncertainty, contrasts, Wald tests, confidence intervals, and exponentiated effects. `GeeDiagnostics` reports cluster scores, leverage, Cook distances, and one-step or exact deletion fits.

Select an adjustment explicitly:

```java
GeeOptions adjusted = GeeOptions.builder()
    .correlation(GeeCorrelation.AR1)
    .covariance(GeeCovariance.BIAS_CORRECTED)
    .method(GeeMethod.BIAS_REDUCED)
    .inference(GeeInference.CLUSTER_T)
    .estimateDispersion()
    .build();
```

Prepare and sort once when fitting several working structures or warm starts:

```java
PreparedGeeData prepared = Gee.prepare(
    response, rowMajorDesign, rows, columns, subjectId, visit,
    weights, offset, preparationOptions, family);

GeeResult fit = Gee.fit(prepared, family,
    preparationOptions.toBuilder()
        .correlation(GeeCorrelation.AR1)
        .initialCoefficients(previousFit.coefficients())
        .parallelism(8)
        .build(),
    BackendPolicy.CPU);
```

Positive-definite projection is enabled by default for unstructured, fixed, and user-designed working matrices. Set `positiveDefiniteProjection(false)` to fail instead. `associationTolerance`, `scaleTolerance`, `scoreTolerance`, and `associationDamping` provide independent convergence control.

## Validation and performance

`GeeRReferenceTest` compares Gaussian, binomial, Poisson, weighted/offset, irregular-wave, and standard working-structure fits with `geepack` 1.3.13 and `geer` 0.1.0. Coefficients and robust covariance use tight or method-appropriate tolerances; nuisance estimates and small-sample covariance use broader tolerances where package conventions differ. The generator records exact package versions and every reference value.

`GeeBenchmark` and `src/benchmark/r/gee_benchmark.R` use the same deterministic design. On the development machine, preparation of 50,000 observations took 0.09 seconds. Eight-thread Java fits took 1.00 seconds for exchangeable, 1.11 seconds for AR(1), and 2.05 seconds for unstructured, versus 0.62, 0.87, and 1.78 seconds for geepack. Java and R coefficients agreed closely. These are environment-specific smoke numbers, not universal performance claims; parallelism helps primarily when there are many clusters and enough per-cluster work.

For very small numbers of clusters, prefer a leverage correction or `JACKKNIFE` with `CLUSTER_T`, and report the number of clusters. Exact deletion performs one additional fit per cluster and is therefore opt-in through covariance selection or `exactClusterDeletion(true)`. A working structure improves efficiency only when it is a useful approximation; robust covariance protects coefficient inference from correlation misspecification but not from a misspecified marginal mean model.

`NominalGee` currently uses working independence for the mean iteration and cluster-robust covariance for dependence. Ordinal local odds ratios operate on the cumulative binary representation. Neither API should be described as a full replication of multgee's joint multinomial local-odds-ratio solver.
