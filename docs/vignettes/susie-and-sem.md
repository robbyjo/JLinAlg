# SuSiE fine mapping, colocalization, and structural equation models

> **Performance status:** SuSiE and multi-signal colocalization are validated
> against susieR and coloc. SuSiE is benchmarked on the official
> `N3finemapping` vignette data, and SEM is validated against `lavaan` and
> benchmarked on a TOPMed cardiometabolic path model.

For a file-based executable-JAR workflow, use the
[CLI-only SuSiE tutorial](cli-susie-tutorial.md), then the
[CLI-only colocalization tutorial](cli-colocalization-tutorial.md). The Java
examples below additionally cover individual-level and sufficient-statistic
SuSiE plus SEM. These CLI commands are included in v0.3.5 and later.

## SuSiE with individual-level data

Rows are samples and columns are candidate variables in the same order as the
names list:

```java
SusieOptions options = new SusieOptions(
    5,      // maximum single effects
    200,    // iterations
    1e-8,   // convergence tolerance
    0.2,    // absolute prior effect variance on standardized X
    true,   // estimate residual variance
    0.95,   // credible-set coverage
    0.5);   // minimum credible-set purity

SusieResult fineMap = Susie.fit(
    phenotype, genotypeMatrix, variantNames,
    options, BackendPolicy.CPU);

if (!fineMap.converged()) {
    throw new IllegalStateException("SuSiE did not converge");
}
double[] pip = fineMap.pip();
double[] posteriorBeta = fineMap.posteriorMean();
List<CredibleSet> sets = fineMap.credibleSets();
```

Individual-level inputs are centered and scaled once. Posterior means are
returned on the original predictor scale with a fitted intercept.

## Summary z scores and LD

```java
double[] z = {12.0, 0.2, -0.1, -10.0, 0.3};
double[][] ld = {
    {1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0},
    {0, 0, 1, 0, 0},
    {0, 0, 0, 1, 0},
    {0, 0, 0, 0, 1}
};

SusieResult summary = Susie.fitSummary(
    z, ld, 1_000, List.of("a", "b", "c", "d", "e"),
    new SusieOptions(2, 200, 1e-8, 0.2, false, 0.95, 0.5),
    BackendPolicy.CPU);
```

Summary mode uses the same finite-sample z transformation as `susieR::susie_rss`.
It requires an integer sample size greater than two, unique variant names,
and a finite symmetric positive-semidefinite LD correlation matrix from an
ancestry-matched reference. Singular LD is allowed. External LD is checked
with scaled pivoted Cholesky: O(P³) worst-case work and O(P²) storage, unlike
susieR's default `check_input = FALSE`. This extra safety cost is included in
the new summary-mode benchmark; it does not apply to the individual-level
cross-product builder, whose matrix is already known to be PSD.
Use `fitSufficientStatistics` when `X'X`, `X'y`, and `y'y` are available on a
different scale. That method checks the joint cross-product matrix for PSD,
so an impossible `X'y` is rejected as well. The summary z transformation uses
`hypot` to avoid squaring overflow, and tiny nonzero PIPs are retained using
log-complement accumulation.

The conjugate update also retains representable posterior means for tiny
positive priors: it forms the posterior variance on the smaller-variance
scale and evaluates products/quotients with separated binary exponents.
For scalar `X'X=1e-200`, `X'y=1e-100`, `y'y=1`, `n=2`, and fixed prior
variance `1e-200`, the posterior mean is `1e-300`, not a zero caused by an
underflowed shrinkage ratio.

An exhausted iteration budget returns `converged() == false` and retains the
residual variance used to produce the returned posterior and objective; it
does not report an extra unperformed variance update. A large objective
decrease cannot satisfy the convergence test. Check convergence before using
the credible sets downstream. Credible-set purity is the minimum absolute within-set LD;
always inspect it alongside coverage and PIP.

On the bundled 574-by-1,001 vignette benchmark, the portable CPU path matches
fixed-prior susieR PIPs and coefficients within `2e-10` and is 4.83 times faster
on the documented i9-9900K run. See
[performance benchmarks](../performance-benchmarks.md) for exact commands and
environment.

The [boundary audit](../../src/benchmark/resources/genetic-audit/AUDIT.md)
preserves the existing N3 reference tolerances and adds exact one-effect
normal-normal posterior comparisons, indefinite/singular LD cases, extreme
z scores, and iteration-limit consistency. Its base-R one-effect timing is a
closed-form reference, not a fresh `susieR` package speed comparison; Java
runs IBSS to convergence and validates the external LD in that workload.

## Colocalize SuSiE signals

Two JLinAlg fits can be colocalized directly:

```java
ColocSusieResult coloc = ColocSusie.analyze(exposureFineMap, outcomeFineMap);

for (ColocSignalPair pair : coloc.signalPairs()) {
    System.out.printf("effects %d/%d: H4=%.4f%n",
        pair.trait1EffectIndex(), pair.trait2EffectIndex(),
        pair.posteriorH4());
}
double[] sharedVariantPosterior = coloc.sharedVariantPosterior(0);
```

`ColocSusieInput` also accepts an L-by-P matrix of log Bayes factors, which is
the direct equivalent of selected rows from `susieR`'s `lbf_variable` output.
Variant IDs are intersected in the first input's order. Defaults match
`coloc.susie`: p1 and p2 are 1e-4, p12 is 5e-6, and signal pairs with less
than 0.5 posterior overlap are omitted. `ColocOptions` can disable trimming,
change these priors, or provide positive per-variant prior weights.

The regression fixture is generated from coloc 5.2.3's bundled
`coloc_test_data` with susieR 0.14.2. Both the single-signal D1/D2 example and
the multi-signal D3/D4 example match R's H0-H4 and conditional-H4 variant
posteriors to numerical precision. Regenerate it with
`src/test/resources/r-reference/generate-coloc-susie-reference.R`.

## Observed-variable SEM from covariance

This path model estimates `x -> m -> y` and free residual variances:

```java
SemModel model = SemModel.builder("x", "m", "y")
    .regression("m", "x", 0.4)
    .regression("y", "m", 0.6)
    .variance("x", 1.0)
    .variance("m", 0.8)
    .variance("y", 0.5)
    .build();

double[] covariance = {
    1.0, 0.5, 0.35,
    0.5, 1.0, 0.70,
    0.35, 0.70, 1.0
};

SemFitResult sem = Sem.fitCovariance(
    covariance, 2_000, model,
    SemOptions.defaults(), BackendPolicy.PREFERRED);

SemParameterEstimate path = sem.parameter("y~m");
System.out.printf("beta=%g se=%g z=%g p=%g%n",
    path.estimate(), path.standardError(),
    path.zStatistic(), path.pValue());
```

Results include the likelihood chi-square, CFI, TLI, RMSEA, SRMR, AIC, BIC,
and implied covariance. A just-identified model can fit perfectly by
construction; fit indices are informative only with positive model degrees of
freedom.

## SEM from rows and equality constraints

```java
SemModel constrained = SemModel.builder("x", "m", "y")
    .regression("equalPath", "m", "x", 0.4)
    .regression("equalPath", "y", "m", 0.4)
    .build();

SemFitResult fromRows = Sem.fit(
    observedRows, constrained,
    SemOptions.defaults(), BackendPolicy.PREFERRED);
```

Repeated labels impose equality constraints. Use the builder's fixed
regression/variance/covariance methods for known parameters. The row-data fit
uses complete-case covariance ML according to `SemOptions`.

This compact example uses observed-variable covariance ML. The separate
joint APIs also support latent measurement loadings and structural means,
Gaussian pattern FIML, ordinal probit pairwise likelihood, case/cluster
sandwich corrections, continuous-model modification indices, and indirect
delta-method inference. See the [current SEM vignette](sem.md) and
[advanced-method validation](../advanced-validation.md) for executable examples
and R comparisons. Available-pair missing ordinal PML under MCAR is now
available through `SemOrdinal.fitPairwiseMissing`. Ordinal PML is not DWLS/WLSMV;
mixed continuous/ordinal models, general MAR ordinal missingness,
multigroup invariance, ordinal modification indices, and robust scaled FIML
fit statistics remain open. The [TOPMed report](../topmed-sem-performance.md)
records the historical observed-variable benchmark, not these newer workloads.
