# Structural equation modeling

`Sem` fits observed and latent RAM models jointly by Gaussian maximum
likelihood. With directed paths `A`, disturbance covariance `S`, structural
intercepts `alpha`, and the selector `F` for observed variables, the model is

```
T = inverse(I - A)
mu = F T alpha
Sigma = F T S T' F'
```

Loadings, latent regressions, disturbance covariances and intercepts enter one
likelihood. No PCA scores are substituted for latent variables. Input columns
follow `model.variables()`, which contains only observed variables;
`model.latentVariables()` lists the additional RAM variables.

## Joint latent measurement and structural paths

```java
SemModel model = SemModel.builder("x1", "x2", "x3", "y1", "y2", "y3")
    .latent("f", "g")
    .fixedLoading("x1", "f", 1)
    .loading("l2", "x2", "f", .8)
    .loading("l3", "x3", "f", 1.1)
    .fixedLoading("y1", "g", 1)
    .loading("l5", "y2", "g", .9)
    .loading("l6", "y3", "g", .7)
    .regression("b", "g", "f", .5)
    .variance("vf", "f", 1)
    .variance("vg", "g", .6)
    .meanStructure()
    .build();

SemFitResult fit = Sem.fit(rows, model);
if (!fit.converged() || !fit.informationAvailable()) {
    throw new IllegalStateException("Check convergence and model identification");
}
double loading = fit.parameter("l2").estimate();
double[] means = fit.impliedMeans();
double[] covariance = fit.impliedCovariance();
double[] parameterCovariance = fit.parameterCovariance();
```

The first loading of each factor fixes its scale here. Alternatively, fix a
factor variance and estimate all loadings. Identification is the caller's
responsibility: an unfixed factor scale, rotational freedom or confounded
paths can make information singular even when the numerical score is small.
Identification is checked on the analytic observed-distribution Jacobian using
pivoted, reorthogonalized QR of normalized parameter columns. Singular or
numerically weak information returns `NaN` covariance entries and
`informationAvailable() == false`; it is never repaired by adding a ridge.
For these fits, inferential degrees of freedom are `-1`, and chi-square,
p-value, CFI, TLI, RMSEA, AIC and BIC are unavailable. In particular, a redundant
latent mean must not change significance merely by adding a parameter label.
`fitTestsAvailable()` distinguishes available global tests from unavailable ones.

Each unspecified observed or latent disturbance variance defaults to a free
parameter starting at one. Free variances use log coordinates. Fixed zero
disturbance variances are allowed when the overall observed covariance is
positive definite. Indefinite disturbance matrices and singular path systems
are rejected. Repeated parameter labels impose equality constraints, including
shared loadings and paths. Labels cannot mix variances with other parameter
types because their optimization coordinates differ.

## Means, intercepts, and sufficient statistics

`meanStructure()` adds free observed intercepts where unspecified. Latent
intercepts default to fixed zero. Use `intercept`, `fixedIntercept`, and their
labelled forms to specify a different identified structure. The reported
intercepts are equation intercepts, and `impliedMeans()` returns `F T alpha`.
For example, in `m = alpha_m + a*x + error`, `alpha_m` generally differs from
the sample mean of `m`.

```java
SemModel mediation = SemModel.builder("x", "m", "y")
    .regression("a", "m", "x", .4)
    .regression("b", "y", "m", .6)
    .regression("direct", "y", "x", .1)
    .meanStructure()
    .build();
SemFitResult fit = Sem.fit(rows, mediation);
SemInference.IndirectEffect indirect = SemInference.indirect(fit, .95, "a", "b");
```

Without an explicit mean structure, complete rows are centered and the model
fits their ML covariance. `fitCovariance` preserves this covariance-only
contract. `fitMoments(covariance, means, n, model, options, backend)` fits supplied
means and a covariance to an explicit mean/intercept model. Covariances are
dense row-major and use divisor `N`, not `N-1`. The full parameter covariance
uses the order of `fit.parameters()` and is transformed back from log variance
coordinates, including all off-diagonal terms. Complete-data standard errors
use expected normal information.

## Constrained FIML with missing observations

```java
SemFimlResult missing = SemFiml.fit(rowsWithNaN, model);
SemFitResult fit = missing.fit();
```

For each missingness pattern, the engine extracts the model-implied observed
mean and covariance and minimizes that pattern's Gaussian likelihood. Pattern
counts, means and centered cross-products are accumulated once; every
objective evaluation uses the appropriate observed submatrix. This is the
constrained model's observed-data likelihood, including the normalizing
constant. It does not estimate a saturated EM covariance and feed that
covariance to complete-data ML.

FIML adds free observed intercepts when the supplied model has no mean
structure. Explicit intercept constraints are preserved. Entirely missing
rows contribute no likelihood, are excluded from `observations()`, and do not
add a missingness pattern. Each variable must have at least two observations;
infinite values are errors. `SemFimlResult.means()` and `covariance()` are the
constrained model's implied moments. Its legacy `iterations()` field now
reports objective evaluations, as does `fit.functionEvaluations()`.

FIML parameter covariance inverts the observed likelihood Hessian, computed
by centered differences of analytic scores. Saturated and independence
models are separately optimized over the same missingness patterns to obtain
the likelihood-ratio chi-square, CFI, TLI and RMSEA. Their fits never determine
the constrained model's parameter estimates. If a reference fit fails or its
information is singular, fit indices remain `NaN`. Auxiliary moments are merged
using centered within-/between-pattern moments, avoiding cancellation at large
locations. An auxiliary failure never discards the target model's valid
likelihood or estimates. This includes complete-data models whose constrained
likelihood is regular even though the centered sample covariance is singular.
SRMR uses H1/sample variances for standardization. For a mean structure it
includes squared standardized mean residuals and divides by
`p*(p+1)/2 + p`; covariance-only models divide by `p*(p+1)/2`.

The usual ignorable missingness assumptions, including MAR with distinct
missingness-model parameters, are required for FIML interpretation.
`MissingDataPolicy.OMIT` on `Sem.fit` continues to select complete rows and
does not invoke FIML.

## Robust inference and modification indices

```java
SemInference.RobustResult robust = SemInference.robust(fit);
SemInference.RobustResult clustered = SemInference.robust(fit, clusterIds);
SemInference.IndirectEffect robustIndirect = SemInference.indirect(
    fit, clustered.parameterCovariance(), .95, "a", "b");

var changes = SemInference.modificationIndices(fit,
    SemInference.Modification.regression("y", "x"),
    SemInference.Modification.covariance("x1", "x2"));
```

Supply only modifications appropriate to the fitted model; the example
illustrates the two candidate constructors. The no-candidate overload scans
omitted observed disturbance covariances. A candidate-specific call can also
release a fixed path or intercept. It rejects a parameter already free.
The score test uses the Schur complement
`I_cc - I_cu inverse(I_uu) I_uc`, and projects nuisance scores out of the
candidate score. Thus it is not the unadjusted scalar `score^2 / I_cc` helper.
The result includes the one-df statistic, expected parameter change and
p-value. A release that is unidentified has `NaN` statistics. Releasing an
equality constraint among already free parameters is not implemented by this
candidate API.

Robust inference derives individual Gaussian likelihood scores from the fit.
Complete ML uses expected-information bread, matching lavaan MLM; FIML uses
observed-information bread. Cluster covariance sums scores within clusters
before forming the meat and uses CR0, with no small-sample multiplier. IDs
must follow the informative rows retained by the fit. Fewer than two clusters
are rejected. Robust covariance requires row data; a covariance matrix alone
does not contain empirical fourth moments.

Complete-data robust results include a Satorra-Bentler mean scaling correction
from the projected empirical moment covariance; clustering uses cluster sums
in the same projection. FIML now uses an observed-information mean-scaling
projection with centered case scores. Unresolved/nonpositive scaling suppresses
the corrected statistic. This is not advertised as lavaan's default MLR/Mplus
Yuan-Bentler, MLMV, or finite-cluster corrections. See the
[new estimator contracts and independent R checks](../estimator-extensions.md).

Product delta inference accepts two or more paths and contracts the entire
parameter covariance with the product gradient. Repeated labels and zero
path estimates are handled without dividing by an estimate. These are
normal/delta intervals, not bootstrap or distribution-of-product intervals.

## Ordinal SEM

`SemOrdinal` jointly estimates probit thresholds and structural parameters by
pairwise maximum likelihood (PML), an
[ordinal estimator also supported by lavaan](https://lavaan.ugent.be/tutorial/cat.html).
It maximizes the sum of bivariate ordinal cell log probabilities, not a
Gaussian fit to category codes or a fit to fixed empirical thresholds.

```java
SemModel ordinalModel = SemModel.builder("z1", "z2", "z3", "z4")
    .latent("f").fixedVariance("f", 1)
    .loading("a", "z1", "f", .8).loading("b", "z2", "f", .8)
    .loading("c", "z3", "f", .8).loading("d", "z4", "f", .8)
    .fixedVariance("z1", 1).fixedVariance("z2", 1)
    .fixedVariance("z3", 1).fixedVariance("z4", 1)
    .build();
SemOrdinal.Result ordered = SemOrdinal.fit(
    categories, new int[] {3, 3, 3, 3}, ordinalModel);
```

Rows contain integer categories `0..K-1`. Observed residual variances must be
fixed positive to identify response scales; the example matches lavaan
`parameterization="theta", std.lv=TRUE, estimator="PML"`. Thresholds are on
that unstandardized latent-response scale. Structural intercepts are excluded
because free thresholds already determine response location. Ordered
threshold gaps use log coordinates. Analytic cell-probability derivatives
propagate through standardized thresholds and model-implied correlations.
Bivariate normal probabilities use the Plackett integral with `rho=sin(t)`
to remove the near-unit-correlation singularity, with adaptive 16/32-point
Gauss-Legendre error checks. Small rectangle probabilities use a positive
conditional-normal integral with stable tail differences, avoiding subtraction
of nearly equal CDFs. Zero-count cells enter neither optimization nor case
scores. Observed composite information and case or cluster
score covariance produce the full Godambe/sandwich parameter covariance,
including threshold-loading cross-covariances. The analytic Jacobian of
standardized thresholds and correlations is checked for rank before inversion;
fixing residual variances alone does not establish identification. Excess
structural parameter counts are rejected. Other rank/information failures set
the entire covariance and every SE/p-value to `NaN`, with
`informationAvailable() == false`, including threshold inference. Sandwiches
are formed from outer products of influence vectors to avoid negative variance
artifacts from dense-matrix cancellation.

`pairwiseLogLikelihood()` is explicitly a composite likelihood, not the
multivariate ordinal full likelihood. The result does not fabricate ML AIC,
global chi-square, or WLSMV corrections. `fit` supports complete all-ordinal
data, including binary indicators and observed ordinal paths.
The separate `SemMixed`, `SemDwls`, and `SemMultigroup` APIs add mixed responses,
joint MAR likelihood, DWLS/WLSMV, ordinal likelihood modification indices and
multigroup invariance. See [estimator extensions](../estimator-extensions.md)
for their dimension limits, moment-covariance schema, and examples. These do
not provide composite likelihood-ratio calibration for PML.
Empty marginal categories must be collapsed explicitly.
Correlations with magnitude at least `.9999` and cell probabilities lost to
floating-point underflow or integration-error failures are rejected rather than
clipped and reported as successful fits. Near-boundary or rare-category problems
may not converge.

### Missing ordinal responses

```java
// Missing categories are -1; observed categories remain 0..categoryCount-1.
var ordinal = SemOrdinal.fitPairwiseMissing(data, categoryCounts, model);
// Optional cluster IDs follow the original input rows.
var clustered = SemOrdinal.fitPairwiseMissing(
    data, categoryCounts, model, SemOptions.defaults(), clusterIds);
```

Each pair contributes only where both responses are observed. Rows with fewer
than two observed responses contribute no pairs and are removed; `observations()`
counts retained rows. Every pair needs at least two jointly observed rows,
and every declared category must occur among the retained observations.
This is available-pair PML: consistency requires MCAR, or a separately justified
pair-specific observation mechanism. It is not an ordinal FIML method for
arbitrary MAR missingness. Threshold and structural scores from the same case
are combined before forming the case/cluster sandwich. Nonconvergence suppresses
the entire covariance and all SEs/p-values in both ordinal fitting APIs.
See [independent validation](../feasible-extensions-validation.md).

## Numerical validation and timing

The R generator and frozen data/results are in
`src/test/resources/r-reference/generate-sem-joint-reference.R` and
`src/test/resources/r-reference/sem-joint/`. They use R 4.6.1, lavaan 0.7.2,
and 600 observations. The continuous fixture fits two measurement factors
with a latent regression and nonnormal disturbances. Additional cases cover
structural intercepts, a constrained intercept, MAR missingness, mediation,
MLM covariance, cluster CR0 covariance and efficient modification indices.
The ordinal fixture uses four three-category indicators with jointly fitted
thresholds. Every entry of parameter covariance is compared after mapping
parameter labels. The PML likelihood is independently evaluated from R's
bivariate normal cell probabilities; lavaan's internally scaled optimizer
objective is not treated as a log likelihood.

`SemJointTest` also checks the complete-row covariance-inflation counterexample,
empty rows, analytic RAM scores against finite differences, singular
information, negative-correlation PCA initialization, and an exact negative
binary probit model. It retains the earlier observed-path lavaan fixture.

`SemReviewRegressionTest` promotes the six independent cross-review defects to
tracked regressions. Its R generator is
`src/test/resources/r-reference/generate-sem-review-reference.R`; frozen checks
are in `sem-review.properties` in the same directory. They include the exact
orthant identity at correlations through `+/- .99989`, its correlation
derivative, rectangle probabilities as small as `7.36e-74` from independent
R conditional-normal integration, mean-inclusive lavaan SRMR, ordinal rank
failure, and regular target likelihoods with problematic auxiliary H1 models.
The isolated JUnit run passes 33 cases, including the separate strict-backend
policy regressions. The Java timings below were refreshed after these repairs;
R timings are the saved five-batch measurements for the unchanged estimands.

The refreshed September 8, 2026 run on this shared Windows host produced the
following timings, after the adaptive ordinal, mean SRMR and identification repairs.
Each median is from five batches of five fits, following warm-up. These are
small-model measurements under concurrent host load, not universal speed
guarantees. Java's FIML timing includes its fitted H1 and independence models.

| Case | Maximum estimate error | Maximum covariance error | Java seconds | R seconds | R / Java |
|---|---:|---:|---:|---:|---:|
| Continuous latent | 2.41e-8 | 2.10e-10 | .00095678 | .026 | 27.2 |
| Latent with means | 2.39e-8 | 1.49e-10 | .00098646 | .026 | 26.4 |
| Constrained FIML | 7.34e-8 | 1.65e-9 | .00649398 | .112 | 17.2 |
| Complete MLM inference | 2.39e-8 | 1.87e-10 | .00233718 | .028 | 12.0 |
| Ordinal PML | 2.63e-8 | 3.77e-9 | .00584022 | .378 | 64.7 |

Absolute likelihood discrepancies were below `9e-12`. The Java benchmark
rejects unconverged fits and accuracy failures before timing, and consumes
estimates, covariance and likelihood through a checksum. R checks convergence
for every timed fit. R's higher-level fitting and output preparation do more
work than these small Java APIs; the estimands and reported likelihoods agree,
but the timings do not imply full feature parity.

Regenerate the R fixtures from the repository root:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' `
  src/test/resources/r-reference/generate-sem-joint-reference.R
```

Build and run with the dependency version pinned by Gradle:

```powershell
.\gradlew.bat benchmarkSemJoint
```

The joint engine currently uses portable Java CPU matrix kernels. CPU,
PREFERRED and AUTO policies are supported; explicit accelerator/native policies
throw `UnsupportedOperationException` before fitting, preserving strict selection
semantics. No accelerator speedup is claimed. Analytic covariance derivatives
and cached missing-pattern sufficient statistics reduce work independently
of the observation count during optimization. Information and final robust
case scores incur additional work. Very large latent systems and many
distinct missingness patterns have not been performance-certified.

The legacy `LatentMeasurement.fit(data, factorCount)` remains a descriptive PCA
extractor, now with a non-null initialization for negative correlations.
`LatentMeasurement.fit(data, model)` runs joint RAM ML. Likewise,
`SemMeanStructure.fit(data)` remains an intercept-only summary;
`SemMeanStructure.fit(data, model)` jointly fits structural intercepts. These
compatibility summaries are not substitutes for the fitted model APIs.

Multigroup invariance, nonlinear constraints, random effects, finite-sample
cluster degrees of freedom, and automatic identification selection are outside
the current implementation. Log-variance coordinates exclude negative
Heywood variance estimates; solutions approaching zero need substantive and
numerical review. `converged()` certifies the configured per-case score
tolerance, not uniqueness, global optimality, model identification, or good
substantive fit.
