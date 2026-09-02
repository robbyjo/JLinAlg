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

`design` includes the intercept when one is wanted. Cluster identifiers need not be sorted; `visit` supplies the within-cluster wave and must be unique inside a cluster. The overload accepting a flat row-major design avoids allocation of a `double[][]`. `GeeFormula.fit` accepts a reusable `CompiledFormula`.

The result exposes coefficients, selected covariance and Wald inference, naive/robust/DF-adjusted/bias-corrected covariance matrices, fitted means, Pearson residuals, working-association and scale estimates, QIC/QICu/CIC/QICC, retained-row mapping, convergence diagnostics, linear-contrast tests, and backend provenance.

## First-release surface

- GLM families and links through `GlmFamily`, with observation weights, offsets, estimated or fixed dispersion, missing-row omission, and CPU/accelerated backend selection.
- Independence, exchangeable, AR(1), m-dependent, Toeplitz, unstructured, fixed, and user-designed working association structures.
- Correlation and binary alternating-odds-ratio association models.
- Model-based and robust sandwich inference plus QIC-family model-selection criteria.
- Rectangular and contiguous row-major entry points, formula integration, immutable defensive-copy results, and reproducible R-reference tests.

## Extended surface

The five follow-on areas are included:

1. Small-sample inference: cluster degrees-of-freedom scaling and a cluster-leverage-adjusted sandwich covariance.
2. Bias reduction: iterative bias-reduced estimating equations, a one-step bias correction, Jeffreys-penalized updates, and one-step/hybrid Jeffreys modes.
3. Nuisance models: observation-level log-dispersion regression, fixed association matrices, and user-defined association design matrices.
4. Ordinal outcomes: `OrdinalGee.fit` implements a proportional-odds cumulative-logit GEE using cumulative binary indicators and a joint sandwich covariance. Its first version intentionally requires working independence.
5. Clinical reporting: `ClinicalGee.marginalMeans` and `ClinicalGee.contrast` provide reference-grid adjusted means, delta-method uncertainty, Wald tests, confidence intervals, and exponentiated link-scale effects.

Select an adjustment explicitly:

```java
GeeOptions adjusted = GeeOptions.builder()
    .correlation(GeeCorrelation.AR1)
    .covariance(GeeCovariance.BIAS_CORRECTED)
    .method(GeeMethod.BIAS_REDUCED)
    .estimateDispersion()
    .build();
```

## Validation and performance

`GeeRReferenceTest` compares an unequal-cluster Gaussian exchangeable fit directly with `geepack` 1.3.13 and `geer` 0.1.0. Coefficients, association, and robust covariance agree with `geer` to numerical precision on that oracle; the leverage-corrected covariance is checked separately with a method-appropriate tolerance because small-sample conventions differ between packages.

`GeeBenchmark` is a reproducible cluster-scaling smoke benchmark. On the development machine, its 10,000-cluster, 50,000-observation CPU run took about 0.94–1.12 seconds after warm-up for exchangeable, AR(1), and unstructured structures (the first, independence run included class/backend warm-up and took 1.94 seconds). These are environment-specific smoke numbers, not a cross-library benchmark.

For very small numbers of clusters, prefer `BIAS_CORRECTED` or `DF_ADJUSTED` over the asymptotic robust covariance and report the number of clusters. A working structure improves efficiency only when it is a useful approximation; robust covariance protects coefficient inference from correlation misspecification but not from a misspecified marginal mean model.
