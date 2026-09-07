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

The GLS path adds a common random-effects variance to the supplied covariance
and retains the generalized Q diagnostic. For clustered meta-regression, use
the cluster labels and moderator matrix directly:

```java
MetaClusterRobustResult robust = MetaClusterRobust.fit(
    studies, clusterLabels, moderators, List.of("dose"),
    MetaAnalysisOptions.randomEffects(), BackendPolicy.PREFERRED);
```

This reports a study-cluster sandwich covariance. It is deliberately explicit:
it does not treat correlated effects as independent or invent a cluster ID.

## Effect sizes and publication-bias diagnostics

Construct additive-scale inputs before pooling:

```java
MetaEffectSize logOr = MetaEffectSizes.logOddsRatio(20, 80, 10, 90);
MetaStudy study = logOr.study("trial-1");
MetaPublicationBiasResult bias = MetaPublicationBias.diagnose(studies);
```

The diagnostic result contains the Egger intercept test and a deterministic
rank-correlation test. They are screening diagnostics, not proof of selective
publication and not substitutes for a prespecified sensitivity analysis.

## Current scope

The original univariate engine still assumes independent study estimates. The
new correlated and robust APIs are intentionally bounded: correlated pooling
is intercept-only, the random-effects extension uses a common variance added
to the supplied covariance, and publication-bias tests are asymptotic
diagnostics. Full multilevel moderator covariance structures remain future
work. See the [numerical contract](../numerical-contract.md) for exact
formulas.
