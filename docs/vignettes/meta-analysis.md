# Meta-analysis and meta-regression

The meta engine expects an effect estimate and positive sampling standard error
for each independent study. Transform ratio measures to a suitable additive
scale, such as log odds ratios, before fitting and back-transform only for
presentation.

## Define the studies

```java
List<MetaStudy> studies = List.of(
    new MetaStudy("study-a", 0.20, 0.10),
    new MetaStudy("study-b", 0.50, 0.20),
    new MetaStudy("study-c", 0.10, 0.15),
    new MetaStudy("study-d", 0.70, 0.25));
```

The third field is an SE, not a variance. Duplicate names are allowed by the
numerical engine but should normally be resolved in data preparation.

## Fixed-effect pooling

```java
MetaAnalysisResult fixed = MetaAnalysis.fit(
    studies, MetaAnalysisOptions.fixedEffect(), BackendPolicy.PREFERRED);

System.out.printf("beta=%g se=%g p=%g Q=%g Qp=%g I2=%g%n",
    fixed.pooledEffectSize(), fixed.standardError(), fixed.pValue(),
    fixed.cochranQ(), fixed.cochranQPValue(), fixed.iSquared());
```

This is inverse-sampling-variance pooling. The confidence interval describes
the common pooled effect. A fixed-effect result intentionally has no
random-effects prediction interval.

## Random-effects pooling

REML is the default tau-squared estimator:

```java
MetaAnalysisOptions randomOptions = MetaAnalysisOptions.builder()
    .method(MetaAnalysisMethod.RANDOM_EFFECT)
    .tauSquaredEstimator(TauSquaredEstimator.REML)
    .inferenceMethod(MetaInferenceMethod.NORMAL)
    .confidenceLevel(0.95)
    .build();

MetaAnalysisResult random = MetaAnalysis.fit(
    studies, randomOptions, BackendPolicy.PREFERRED);

System.out.printf("beta=%g se=%g tau2=%g PI=[%g,%g]%n",
    random.pooledEffectSize(), random.standardError(),
    random.tauSquared(), random.predictionLower(),
    random.predictionUpper());
```

Select `DERSIMONIAN_LAIRD` for the generalized method-of-moments estimator or
`PAULE_MANDEL` for the generalized-Q solution. Show a sensitivity table when
the number of studies is small or conclusions depend on the estimator.

## Knapp-Hartung inference

```java
MetaAnalysisOptions hkOptions = MetaAnalysisOptions.builder()
    .tauSquaredEstimator(TauSquaredEstimator.REML)
    .inferenceMethod(MetaInferenceMethod.HARTUNG_KNAPP)
    .build();

MetaAnalysisResult hk = MetaAnalysis.fit(
    studies, hkOptions, BackendPolicy.PREFERRED);
```

Knapp-Hartung scales the coefficient covariance by residual dispersion and
uses `k-p` Student-t degrees of freedom. The unmodified scale can occasionally
reduce the SE. Use `MODIFIED_HARTUNG_KNAPP` to cap the scale below at one.
`STUDENT_T` uses t inference without the covariance rescaling.

All results retain ordinary, `log10(p)`, and `-log10(p)`:

```java
double p = hk.pValue();
double log10P = hk.log10PValue();
double minusLog10P = hk.negativeLog10PValue();
```

## Meta-regression

Moderator rows must match study order. Do not add an intercept column when
`includeIntercept` is true:

```java
double[][] moderators = {
    {-1.0, 0.0},
    { 0.0, 1.0},
    { 1.0, 0.0},
    { 2.0, 1.0}
};

MetaRegressionResult regression = MetaRegression.fit(
    studies, moderators, List.of("dose", "highRisk"), true,
    randomOptions, BackendPolicy.PREFERRED);

for (int j = 0; j < regression.coefficientNames().size(); j++) {
    System.out.printf("%s beta=%g se=%g p=%g%n",
        regression.coefficientNames().get(j),
        regression.beta()[j], regression.standardErrors()[j],
        regression.pValues()[j]);
}

System.out.printf("QE=%g QEp=%g QM=%g QMp=%g R2=%g%n",
    regression.residualQ(), regression.residualQPValue(),
    regression.moderatorQ(), regression.moderatorQPValue(),
    regression.heterogeneityRSquared());
```

`Q_E` tests residual heterogeneity with fixed sampling-variance weights.
`Q_M` is the asymptotic Wald chi-square test of all moderator coefficients,
excluding the intercept. Heterogeneity R-squared is the nonnegative reduction
in tau-squared relative to the intercept-only model; it is not ordinary OLS
variance explained.

## Correlated effects and robust inference

For repeated outcomes, multiple estimates from one study, or a multilevel
effect structure, provide the study-level sampling covariance explicitly:

```java
MetaCorrelatedResult correlated = MetaCorrelatedAnalysis.fit(
    studies, samplingCovariance, MetaAnalysisOptions.randomEffects(),
    BackendPolicy.PREFERRED);
```

The GLS path estimates `tau2` in `V + tau2 I`. `V` is the authoritative sampling
covariance; its diagonal is not reconstructed from the study SEs. REML minimizes
the restricted Gaussian likelihood, Paule–Mandel solves the generalized residual
Q equation, and generalized DerSimonian–Laird uses `(Q(0)-(k-p))/trace(P(0))`,
truncated at zero. Reported generalized Q uses the sampling covariance at zero
heterogeneity. Normal, residual t, Hartung–Knapp and modified Hartung–Knapp
inference all apply their stated scaling and confidence level.

For clustered meta-regression with correlated sampling errors:

```java
MetaClusterRobustResult robust = MetaClusterRobust.fit(
    studies, clusterLabels, moderators, List.of("dose"),
    samplingCovariance,
    MetaAnalysisOptions.builder()
        .inferenceMethod(MetaInferenceMethod.STUDENT_T).build(),
    MetaRobustCorrection.CR2, BackendPolicy.CPU);
```

The covariance must be block diagonal between independent clusters; dependence
within a cluster is allowed. At least two clusters are required. CR0 is the
unadjusted sandwich. The old overload retains CR0 as its default. Corrections
follow `clubSandwich` naming:

| Correction | Covariance adjustment | Student-t degrees of freedom |
|---|---|---|
| CR0 | None | G-p |
| CR1 | G/(G-1) | G-p |
| CR1P | G/(G-p), matching `metafor::robust(adjust=TRUE)` | G-p |
| CR1S | G(N-1)/((G-1)(N-p)) | G-p |
| CR2 | Symmetric Bell–McCaffrey bias-reduced adjustment | Coefficient-specific Satterthwaite |

CR2 uses the fitted total covariance as its working target. With residual block
covariance `D_i = V_i - X_i B X_i'`, its adjustment is
`A_i = V_i^(1/2) [V_i^(1/2) D_i V_i^(1/2)]^(-1/2) V_i^(1/2)`.
The Satterthwaite calculation contracts the adjusted cluster scores with the
full residual projection; it is not a substitution of `G-1` for every
coefficient. Singular CR2 residual blocks, nonpositive residual cluster DF,
and unsupported Hartung–Knapp-plus-sandwich combinations are rejected explicitly.
Normal inference remains available, and confidence intervals honor the selected
confidence level. One cluster never produces a zero covariance as a valid fit.

## Hierarchical random moderators

The original six-argument `MetaMultilevelRegression.fit` is known-total-covariance
GLS. Its Gaussian ML likelihood includes `N log(2 pi)`. To estimate hierarchical
covariance, use the overload with random-effect terms:

```java
// Rows explicitly contain the random intercept and slope design.
double[][] randomDesign = new double[studies.size()][2];
for (int i = 0; i < studies.size(); i++) {
    randomDesign[i][0] = 1;
    randomDesign[i][1] = moderators[i][0];
}
List<MetaRandomEffect> terms = List.of(
    new MetaRandomEffect("study", clusterLabels, randomDesign,
        MetaRandomEffect.Structure.UNSTRUCTURED));

MetaHierarchicalResult hierarchical = MetaMultilevelRegression.fit(
    studies, moderators, List.of("dose"), samplingCovariance, true, terms,
    MetaMultilevelRegression.Estimation.REML,
    MetaAnalysisOptions.randomEffects(), BackendPolicy.CPU);
if (!hierarchical.converged()) {
    throw new IllegalStateException("Covariance optimization did not converge");
}
double[] studyCovariance = hierarchical.randomCovariances().get(0);
// Row-major: intercept variance, covariance, covariance, slope variance.

MetaClusterRobustResult hierarchicalRobust = MetaClusterRobust.fit(
    studies, clusterLabels, moderators, List.of("dose"), hierarchical,
    MetaAnalysisOptions.builder()
        .inferenceMethod(MetaInferenceMethod.STUDENT_T).build(),
    MetaRobustCorrection.CR2, BackendPolicy.CPU);
```

The marginal covariance is `V + sum Z_g G_g Z_g'`. Each term shares an estimated
`G_g` across its groups. `DIAGONAL` estimates independent random coefficients;
`UNSTRUCTURED` also estimates their covariances. Multiple terms are optimized
jointly. For nested intercepts, supply
`MetaRandomEffect.intercept("study", studyIds)` and
`MetaRandomEffect.intercept("outcome", studyOutcomeIds)`; nested labels must
include their parent label. Distinct terms may also represent crossed effects,
but a robust cluster must encompass all dependencies.

ML and REML are separate estimands. The REML log likelihood uses orthonormal
error contrasts, including the `log|X'X|` constant, matching `metafor::rma.mv`.
Sampling covariance is known and is not multiplied by a free residual scale.
Cholesky coordinates keep every random covariance positive semidefinite,
including boundary solutions. The optimizer uses two deterministic starts,
scaled coordinates, simplex-spread and likelihood checks, and a numerical
stationarity check. A budget exhaustion returns `converged=false`. Linearly
dependent covariance bases are rejected after projection out of the REML fixed
space. This avoids redundant variance components with arbitrary estimates.

`maximumIterations` is the per-coordinate budget for each optimizer start.
The hierarchical overload requires random-effects options with the default
`tauSquaredEstimator=REML`; its separate `Estimation` argument chooses joint
ML or REML. Scalar DL/PM options are rejected because they do not define this
multicomponent covariance fit. Coefficient inference uses plug-in covariance;
the robust overload conditions on the converged fitted covariance and uses
the options' inference method and confidence level.

## Effect sizes and publication-bias diagnostics

Construct additive-scale inputs before pooling:

```java
MetaEffectSize logOr = MetaEffectSizes.logOddsRatio(20, 80, 10, 90);
MetaStudy study = logOr.study("trial-1");
MetaPublicationBiasResult bias = MetaPublicationBias.diagnose(studies);
```

Egger's standardized-effect regression estimates multiplicative dispersion and
uses Student t with `k-2` degrees of freedom. At least three studies and varying
SEs are required. Its intercept equals the SE slope in PET, not PET's intercept.
The rank diagnostic is approximate Spearman correlation with tie-averaged ranks;
it is not the Kendall/Begg statistic returned by `metafor::ranktest`.

`MetaBiasCorrections.pet(studies, .95)` fits `effect ~ SE`, weighted by `1/SE^2`.
`peese` substitutes sampling variance as the predictor. Both expose the two
coefficients, full covariance, residual t tests, and confidence intervals.
Residual dispersion is estimated without a lower bound of one, matching R's
weighted `lm`; these are multiplicative-dispersion regressions, not additive
random-effects regressions or an automatic PET-to-PEESE selection procedure.

`MetaBiasCorrections.trimAndFill` implements Duval–Tweedie L0, R0 and Q0. It
iteratively trims ranked effects, refits the center with the selected pooling
estimator, estimates the number missing from signed absolute ranks, and fills
by reflection of the selected studies with their original sampling variances.
The augmented dataset is then pooled with the selected estimator. Defaults are
REML, L0, and side selection by a meta-regression on SE. The overload accepts an
explicit `Side.LEFT` or `Side.RIGHT`, especially useful when SEs do not vary.
The result exposes imputed count, side, iteration count, augmented studies and
the full adjusted fit. Oscillation, too few retained studies, or an undefined
Q0 square root produces an exception, not a fabricated correction.

These are sensitivity analyses, not proof of selective publication. Hedges g
uses the documented `J ≈ 1-3/(4df-1)` correction and a delta-method sampling
variance; it is not claimed to reproduce every `escalc` variance convention.
Both sample-size addition and reciprocal calculations use floating point,
including for sample sizes above the integer-product overflow boundary.

## Reference accuracy and measured speed

The frozen inputs and 236 reference scalars are in `src/test/resources/meta`.
`generate-metafor.R` regenerates them with R 4.6.1, metafor 5.0-1 and
clubSandwich 0.7.0. It prefixes the local R library path and requires no network
access once dependencies are installed. Tests cover ML/REML likelihoods,
unstructured and diagonal random moderators, nested variance components,
correlated REML/DL/PM, all sandwich corrections, finite DF, PET/PEESE, Egger,
positive-heterogeneity trim-and-fill, and failures at invalid boundaries.

For trim-and-fill the generator refits the identical augmented R dataset with
tight controls: `trimfill` otherwise drops the original fit's controls in its
final refit. It also sets both Fisher-scoring `threshold` and root-finding `tol`
to `1e-12`, so optimizer defaults do not masquerade as estimator differences.
Generalized correlated DL/PM references use explicit R GLS estimating equations;
the corresponding REML reference is `rma.mv`.

On 2026-09-08, 23 isolated meta-package tests passed with Java 17 source checks
(`javac --release 17 -Xlint:all -Werror`) on Java 25. For the 36-effect,
12-cluster fixture, likelihood agreement was within `1e-8`, benchmarked random
covariance error below `7e-7`, and the known-covariance CR2 benchmark agreed within `1e-12`.
These are fixture results, not a bound for arbitrary designs.

| Workflow | Java CPU ms/fit | R ms/fit |
|---|---:|---:|
| Known covariance GLS | 0.114 | 2.60 |
| Random-moderator REML | 4.44 | 37.3 |
| Nested-intercept REML | 2.32 | 11.9 |
| Correlated GLS + CR2 + t inference | 0.210 | 8.50 |
| Fixed-effect L0 trim-and-fill, 12 studies | 0.135 | 11.1 |
| PET, 12 studies | 0.0131 | 0.400 |

Each timing uses ten warmups and 100 measured fits, excludes fixture loading,
checks reference accuracy before timing, and consumes every returned result
in a checksum. Both implementations fit the same likelihoods/estimands, but
their public APIs have different ancillary diagnostics and overhead. These
small dense workflow timings are not large-study scaling or isolated kernel
speed claims. The recorded values and checksums are in
`src/benchmark/resources/meta/parity-timings-2026-09-08.csv`.

Reproduce from the repository root, after compiling the benchmark source:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/meta/generate-metafor.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/meta_parity_benchmark.R 100
java -Xms256m -Xmx1g -cp 'build/classes/java/benchmark;build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar' org.jlinalg.benchmark.MetaParityBenchmark 100
```

The host may require an elevated R invocation to read dependencies in its user
library. For a shared workspace, compile the meta package and this benchmark
into `build/meta-classes`, then put that directory first on the classpath.

## Current scope

The univariate engine assumes independent study estimates. Scalar correlated
pooling is intercept-only; the hierarchical API supplies fixed and random
moderators. Hierarchical estimation is dense, with quadratic covariance storage
and cubic factorization cost. Its multistart local optimizer is not a proof of
a global optimum. Covariance-component profile intervals, variance-parameter
uncertainty corrections to model-based coefficient inference, and the complete
set of spatial/structured covariance families in `rma.mv` are not implemented.
Model-based residual Hartung–Knapp scaling should not be confused with CR2 or
Kenward–Roger. The hierarchical robust convenience overload currently requires
a fixed intercept and matching moderator rows. CR2 rejects singular or
numerically ill-conditioned residual cluster blocks rather than reporting
unstable finite-sample inference.
