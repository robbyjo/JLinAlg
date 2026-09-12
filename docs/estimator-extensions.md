# SEM, time-series, nonparametric, and integration extensions

Validated September 12, 2026. These Java APIs extend the four estimator items
previously listed in TODO. The implementations retain explicit dimension,
identification, convergence, and inferential limits; they do not imply full
parity with every option in the reference packages.

## SEM

`SemMixed.fit(data, categoryCounts, model)` maximizes a joint observed-data
Gaussian/ordinal-probit likelihood. A category count of zero marks a continuous
response; ordinal responses use integer codes `0..K-1`. Missing responses use
`Double.NaN`. All-missing rows are removed, while singleton observations remain
in the likelihood. Conditional Gaussian densities and ordinal rectangle
probabilities marginalize each missingness pattern. Ignorable MAR and distinct
missingness parameters are assumptions; available-pair PML is not substituted.

```java
var model = SemModel.builder("x", "y").meanStructure()
    .intercept("mu", "x", 0).variance("vx", "x", 1)
    .fixedIntercept("y", 0).fixedVariance("y", 1)
    .regression("b", "y", "x", .3).build();
var mixed = SemMixed.fit(data, new int[]{0, 3}, model);
var covariance = mixed.parameterCovariance();
```

Fix ordinal residual variances to positive values and ordinal intercepts to
identify their scales and locations. Thresholds are estimated jointly; reported
thresholds are actual thresholds, with their full transformed covariance.
Observed information uses numerical derivatives. Failed optimization or
unresolved information suppresses covariance. Every marginal category must be
observed. Deterministic rectangle refinement supports at most **four observed
ordinal responses per row**; larger rows fail explicitly. Duplicate response
patterns are compressed for likelihood evaluations. Dense integration remains
expensive near extreme correlations and can fail its refinement check.

`SemMixed.modificationIndices(fit, candidates...)` tests omitted/fixed paths or
covariances using efficient likelihood scores, projecting out all fitted paths
and thresholds. These are joint-likelihood modification indices; they are not
PML or DWLS score-test corrections. They require a converged identified fit.

`SemDwls.fit(moments, model)` estimates by DWLS and returns a full sandwich
covariance plus a scaled-and-shifted WLSMV statistic. `SemDwls.Moments` contains
sample size, category counts, the statistic vector, and **full Gamma**, the
asymptotic covariance of the sample statistics. A diagonal Gamma loses the
cross-covariances necessary for the robust correction. Use the conventional
`n-1` finite-sample scaling, as in the checked lavaan fixture.

Statistics are ordered as each variable's continuous mean or ordinal
standardized thresholds, then the row-wise lower covariance triangle, omitting
ordinal diagonals. Ordinal latent responses have unit variance in these sample
statistics. Cross-covariances with continuous responses remain in continuous
response units. Gamma must follow the same order and be positive definite.

For residual vector `e`, Jacobian `D`, diagonal weights `W`, and `df=m-k`:

```text
T = (n-1) e' W e
U = W - W D (D' W D)^-1 D' W
a = sqrt(tr((U Gamma)^2) / df)
b = df - tr(U Gamma) / a
T_WLSMV = T/a + b
```

`SemDwls.prepare(data, categories, variableNames...)` supplies MAR-capable
statistics and Gamma from a saturated joint mixed likelihood and empirical
case-score sandwich. This preparation differs from lavaan's default two-stage
ordinal sample statistics. For exact comparison with that estimator, import
its sample statistics and Gamma. Preparation inherits the joint-likelihood
dimension and identification limits.

`SemMultigroup.fit(groups, categories, models, sharedLabels, options)` fits
independent groups jointly. An empty equality set is configural invariance;
share loading labels for metric invariance and continuous intercept/ordinal
threshold labels for scalar invariance. Share only selected labels for partial
invariance. Models must define comparable factor scales and identification.
Threshold labels use `variable|t1`, `variable|t2`, etc.; equality constrains
actual thresholds, including when earlier thresholds remain group-specific.
`SemMultigroup.compare` checks identical data/models and nested equality sets
before reporting a conventional likelihood-ratio test. The correct joint
distribution and regular interior identification are required. This is not a
multigroup WLSMV implementation or a scaled-WLSMV difference test.

`SemInference.robust(fiml.fit())` now supplies an observed-information
mean-scaled FIML likelihood-ratio statistic. It uses the structured saturated
moment Hessian, the fitted model Hessian, and centered case scores; cluster
scores are summed when IDs are supplied. The scaling is
`tr((A1^-1 - D A0^-1 D') B1)/df`. Nonpositive or unresolved scaling produces
`NaN` fit statistics without replacing the covariance. The reference uses an
independently differentiated R likelihood. This particular convention is not
advertised as lavaan's default MLR/Mplus Yuan–Bentler shortcut. Robust CFI/TLI,
finite-cluster corrections, and robust nested-model difference tests are not
provided by this result.

## Time series

```java
var regression = ArimaRegression.fit(y, design, new ArimaOrder(1, 1, 1),
    ArimaOptions.defaults());
var forecast = regression.forecast(futureDesign, .95);
var historicalErrors = regression.smoothErrors();
```

The design explicitly includes intercepts, trends, seasonal indicators, or
other regressors; `includeMean`/`includeDrift` options do not add columns.
Regression coefficients are profiled by pivoted QR on innovations from the
same exact ARIMA filter, with the response's missingness mask applied to every
design column. Integrated/seasonal dynamics use exact diffuse initialization.
Rank-deficient transformed designs, including an intercept eliminated by
differencing, are rejected. Design values must be finite at all dates.

The returned GLS coefficient covariance conditions on fitted dynamics.
Forecast bands and smoothed error states also condition on fitted parameters;
they do not incorporate parameter-estimation uncertainty. Future regressors
must be supplied explicitly. Optimizer failure leaves coefficient covariance
unavailable.

`ArimaSmoothing.smooth(y, ar, ma, differencePolynomial, innovationVariance)`
returns historical state means/covariances and signal means/variances. It
conditions on all observed dates, including future observations relative to
the historical date, and handles arbitrary missing dates and diffuse initial
levels. Dense Gaussian conditioning is bounded to **1024 dates and 128 states**;
it uses cubic observation-matrix work. Existing long-series filtering and
forecasting retain their state-sized storage. A scalable exact diffuse backward
smoother is a separate extension.

## Nonparametric inference

`ProductKernelRegression.predict` supports multiple continuous, unordered, and
ordered predictors, with Gaussian, Aitchison–Aitken, and geometric product
kernels. Local linear terms apply to continuous coordinates. Categories must
occur in training data; ordered codes must be integers. Rank-deficient local
designs and queries with no support fail.

`ProductKernelRegression.infer(x, y, query, types, level)` exactly matches the
query's categorical stratum and automatically uses continuous bandwidths
`robustScale * n_stratum^(-1/(d+2))`. This predictor-only rule undersmooths a
second-order local linear estimator: `n h^d` diverges and `n h^(d+4)` vanishes.
It reports HC3 pointwise errors using original-response influences. Validity is
asymptotic under independent observations, smooth means/variances, an interior
query, positive design density, fixed dimension, and nonvanishing category
probabilities. Small local support is rejected. These are not simultaneous
bands or a finite-sample coverage guarantee.

Quantile inference methods remain distinct:

- `fitExactIidKernel(y, x, tau)` automatically selects a normal-reference
  residual-density bandwidth `1.06 * residualSD * n^(-1/5)`. It requires iid
  continuous errors, positive smooth density, finite residual variance, and a
  regular fixed-dimensional design. The explicit-bandwidth overload remains.
- `fitExactCluster` uses supplied consistent conditional densities and CR0
  cluster score covariance. Clusters must be independent and sufficiently
  numerous; arbitrary dependence is allowed within them.
- `fitExactHac` uses supplied densities and Bartlett HAC covariance with a
  caller-selected lag on time-ordered rows. Weak-dependence/mixing conditions
  and a suitably increasing lag are required for general serial dependence.
- `marginalInterval` uses exact binomial order-statistic inversion for iid
  marginal quantiles, conservatively valid at response mass points. Infinite
  endpoints report insufficient tail information. This does **not** supply
  discrete-response quantile-regression coefficient inference.

## Non-tensor Gaussian random-effect integration

Pass `MultidimensionalQuadratureOptions.nonTensorDefaults()` to the existing
`MultidimensionalGlmmQuadrature.evaluate` or `fit` API. The option selects eight
randomly shifted Halton replicates, with fixed seeds/common random numbers.
A defensive mixture of the mode-centered Gaussian proposal and the Gaussian
prior keeps importance weights bounded for the supported binomial-logit and
Poisson-log likelihoods. Correlated prior precision bases remain supported.

Here `initialOrder`/`maximumOrder` mean points per replicate, and reported nodes
are eight times that count. Convergence requires two consecutive refinements
to satisfy both the change in log likelihood and a three-standard-error
replicate estimate. This is an **estimated integration error**, not a rigorous
error bound or an optimization/global-maximum certificate. Inadequate budgets
return nonconvergence. Dense random-effect operations are bounded to 256
coefficients; numerical integration can fail much earlier. Tensor Hermite
remains the default and retains its existing node budget. Neither method
replaces sparse methods for large structures.

## Validation and reproduction

Checked-in fixtures and their generator are under
`src/test/resources/estimator-extensions` and
`src/benchmark/r/estimator_extensions_reference.R`. R uses the optional
repository-local `build/r-library`; regeneration needs lavaan and numDeriv.
Ordinary Java tests need neither R nor network access once dependencies exist.

The full `check javadoc executableJar` gate passed: **758 tests discovered,
755 passed, three optional CHOLMOD skips, no failures or errors**. Javadoc,
website checks, and the executable JAR also completed. The new extension tests
contribute 21 cases to that suite.

Validation includes lavaan loadings, WLSMV adjustment and full loading
covariance; independently written mixed MAR likelihood and covariance;
observed-information FIML scaling; analytic trivariate/four-dimensional
orthants; actual-threshold partial invariance; R ARIMA regression and historical
states; analytic random-walk/AR bridges; independent R local-linear HC3
algebra; exact enumeration of discrete-quantile coverage; and independent
scalar/eight-dimensional integration identities and exhausted-budget tests.

Reference-gated benchmarks on this Windows host, three warmups and ten timed
runs, executed without a concurrent test workload:

| Workload | Mean time | Checked quantity error |
|---|---:|---:|
| Mixed MAR likelihood, 140 rows | 7.325 ms | log likelihood `1.40e-11` |
| WLSMV, 14 sample statistics | 7.226 ms | adjusted statistic `1.31e-9` |
| AR(1) regression, 80 dates | 2.702 ms | log likelihood `1.73e-12` |
| Diffuse historical smoothing, 80 dates | 1.956 ms | first state `2.54e-6` versus R's finite diffuse initialization |
| Defensive non-tensor integration, 2 dimensions | 369.478 ms | log likelihood `3.04e-6` versus Hermite |
| Mixed kernel prediction, 120 rows | 0.145 ms | linear reproduction `2.66e-15` |

These small workloads are correctness-gated measurements, not general speed or
coverage claims. Non-tensor integration is not expected to beat Hermite at two
dimensions. Reproduce with:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/estimator_extensions_reference.R
.\gradlew.bat check javadoc executableJar --console=plain
.\gradlew.bat benchmarkEstimatorExtensions --console=plain
```

Method references: [lavaan categorical estimators](https://lavaan.ugent.be/tutorial/cat.html),
[multigroup constraints](https://lavaan.ugent.be/tutorial/groups.html),
[Yuan–Bentler implementation](https://github.com/yrosseel/lavaan/blob/master/R/lav_test_yuan_bentler.R),
[R ARIMA regression](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/arima.html),
[R kernel bandwidths](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/bandwidth.html),
and [mixed product kernels](https://stat.ethz.ch/CRAN/web/packages/np/refman/np.html).
