# JLinAlg

**[Explore the project website](https://robbyjo.github.io/JLinAlg/)** for the
complete feature map and browser-friendly worked vignettes.

JLinAlg implements Java linear and mixed-model algorithms on top of
[JDistlib 0.10.0](https://github.com/robbyjo/JDistlib/releases/tag/v0.10.0).
The current release slice provides ordinary least squares (OLS), generalized
linear models (GLMs), dense Gaussian restricted maximum likelihood (REML),
pedigree animal-model REML, and penalized-quasi-likelihood generalized linear
mixed models (GLMM PQL). Gaussian ridge, LASSO, and elastic-net regression and
a first summary-statistics Mendelian-randomization (MR) layer are also included.

## Requirements and build

- A JDK 17 or newer.
- Network access on the first build. The build downloads the pinned
  `jdistlib-all-0.10.0.jar` and verifies its SHA-256 digest before compiling.

On Windows:

```powershell
.\gradlew.bat check
```

On Linux or macOS:

```shell
./gradlew check
```

No system Gradle installation is needed.

## Vignettes and worked examples

The [vignette index](docs/vignettes/README.md) provides end-to-end examples for
every implemented feature group: OLS/GLM and penalized regression, REML/LMM,
pedigree and GLMM PQL, formulas/backends, association and GWAS/TWAS, Mendelian
randomization, meta-analysis/meta-regression, time series, SuSiE, and SEM. Each
vignette explains input layout, result interpretation, performance choices,
and estimator limitations rather than presenting code without its statistical
contract.

The same material is organized as a responsive, system-themed website under
[`site/`](site/index.html). GitHub Pages deploys that directory after changes
land on `main`; visitors can follow their operating-system theme or select an
explicit light or dark preference.

## OLS

```java
double[] y = {1, 2, 5, 7, 9};
double[][] x = {
    {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}
};

OlsResult fit = Ols.fit(y, x);
System.out.println(Arrays.toString(fit.coefficients()));
```

Full-column-rank models use JDistlib's column-pivoted QR decomposition.
Rank-deficient designs fail by default. Set
`RankDeficiencyStrategy.MINIMUM_NORM` to request an SVD-based Moore-Penrose
solution explicitly. Coefficients and individual tests from that solution
depend on the minimum-norm identifying convention and must not be interpreted
as uniquely estimable original parameters.

`OlsResult` includes coefficients, fitted values, residuals, residual variance,
the coefficient covariance matrix, standard errors, t statistics, two-sided
p-values, confidence intervals, numerical rank, likelihood, and backend
provenance.

All fitted-model result types expose a common `associationStatistics()` view.
It returns beta, standard errors, Wald statistics, denominator degrees of
freedom, two-sided p-values, and the statistic's reference distribution. OLS
uses the conventional `N - rank(X)` residual DF, while non-Gaussian standalone
GLMs use asymptotic Wald z inference.

## GLM

Gaussian/identity, binomial/logit, Poisson/log, Gamma/log,
inverse-Gaussian/log, fixed-size negative-binomial/log, quasi-binomial, and
quasi-Poisson GLMs are available through a shared IRLS engine:

```java
GlmResult fit = Glm.fit(y, x, GlmFamilies.binomial());
```

The full overload accepts positive prior weights and an additive linear-
predictor offset. Binomial responses are proportions in `[0,1]`; prior weights
represent trial counts. Poisson responses are nonnegative integer counts, and
an exposure can be represented by its logarithm as the offset.

OLS has the same positive-weight and additive-offset inputs. Both OLS and GLM
default to `MissingDataPolicy.ERROR`; select `OMIT` to compact complete cases
once before numerical fitting. Results report retained original row indices.
`testContrast()` provides joint finite-DF F tests for OLS and asymptotic
chi-square Wald tests for GLMs.

IRLS uses weighted pivoted QR and deviance-reducing step halving. The result
contains fitted means, linear predictors, coefficient covariance, Wald
statistics, confidence intervals, deviance and Pearson residuals, dispersion,
log likelihood/AIC when defined, convergence status, and backend provenance.

## Ridge, LASSO, and elastic net

`PenalizedRegression` fits Gaussian penalized models by cyclic coordinate
descent. The intercept is not penalized; predictors are standardized during
optimization by default and coefficients are returned on their original scale.
Positive observation weights and per-coefficient penalty factors are supported.

```java
PenalizedRegressionResult ridge =
    PenalizedRegression.ridge(y, x, 0.1);
PenalizedRegressionResult lasso =
    PenalizedRegression.lasso(y, x, 0.1);

ElasticNetOptions options = ElasticNetOptions.builder()
    .alpha(0.5)
    .build();
PenalizedRegressionResult elasticNet =
    PenalizedRegression.fit(y, x, 0.1, options);
```

`automaticPath` creates a descending log-spaced lambda path and uses warm
starts. `PenalizedRegressionCrossValidation` performs reproducible K-fold
selection and returns both the minimum-error lambda and the conservative
one-standard-error choice. Preprocessing is learned separately inside each
training fold.

Ridge coefficient beta, SE, effective-residual DF, t statistic, and p-value
are available from `PenalizedRegressionInference.ridge`. For LASSO and elastic
net, `refitActiveSet` provides an optional OLS refit and association results
conditional on the selected active set. Those p-values are not adjusted for
variable selection and should not be presented as selection-valid inference.

## REML

REML fits the Gaussian variance-component model

```text
y ~ N(X beta, V)
V = variance[0] K[0] + ... + variance[q-1] K[q-1]
```

Every `K` is a caller-supplied symmetric positive-semidefinite covariance
basis matrix. Include an identity component for independent residual variance:

```java
List<VarianceComponent> components = List.of(
    new VarianceComponent("group", observations, groupRelationship),
    VarianceComponent.identity("residual", observations)
);

RemlResult fit = Reml.fit(y, x, components);
```

The implementation maximizes the exact dense restricted Gaussian likelihood.
Set `VarianceEstimation.ML` in `RemlOptions` for profile ML and fixed-effect
model comparison.
Variance parameters are optimized on the log scale using Fisher scoring,
step limiting, positive-definite regularization, and likelihood-increasing
backtracking. The result reports convergence explicitly; callers must check
`converged()` before scientific interpretation.

The dense covariance-component API is intentionally general enough for random
intercepts, repeated measures, kinship matrices, and other covariance bases.

Common uncorrelated random-effect terms can be constructed without manually
forming `K`:

```java
VarianceComponent intercept =
    VarianceComponent.randomIntercept("site", siteIds);
VarianceComponent slope =
    VarianceComponent.randomSlope("site:age", siteIds, age);
```

REML fixed-effect association tests default to the requested fast denominator
DF approximation `N - rank(X) - 1`. Opt into coefficient-specific
Satterthwaite DF with:

```java
RemlOptions options = RemlOptions.builder()
    .degreesOfFreedomMethod(DegreesOfFreedomMethod.SATTERTHWAITE)
    .build();
```

The Satterthwaite calculation uses the REML variance-parameter information and
analytic derivatives of `vcov(beta)`. It does not refit the model per
coefficient, but it is more expensive than the default approximation.

Kenward-Roger small-sample inference is also available:

```java
RemlOptions options = RemlOptions.builder()
    .degreesOfFreedomMethod(DegreesOfFreedomMethod.KENWARD_ROGER)
    .build();
```

This option applies the Kenward-Roger bias adjustment to the fixed-effect
covariance before calculating SE, t statistics, and p-values, and reports
coefficient-specific denominator DF. `fixedEffectCovariance()` remains the
unadjusted model covariance; `fixedEffectInferenceCovariance()` returns the
matrix actually used for inference. Kenward-Roger is substantially slower than
the default approximation because it uses dense observation-scale derivatives.

## General linear mixed models

`LinearMixedModel` adds a term-oriented layer over exact Gaussian REML. It
supports multiple independent random-effect terms, including crossed/nested
random intercepts and independent random slopes:

```java
List<RandomEffectTerm> random = List.of(
    RandomEffectTerm.randomIntercept("site", siteIds),
    RandomEffectTerm.randomSlope("site:age", siteIds, age),
    RandomEffectTerm.randomIntercept("batch", batchIds)
);

LinearMixedModelResult fit = LinearMixedModel.fit(y, x, random);
RandomEffectEstimates siteEffects = fit.randomEffects("site");
```

The grouped random-effect design is stored internally in CSR form (one nonzero
per observation for random intercepts or slopes), while `design()` remains a
dense compatibility view. The result contains fixed-effect association
inference, variance components, conditional random-effect modes,
prediction-error variances, conditional fitted values, and conditional
residuals. For independent residual errors, conditional modes and PEVs are
computed through Henderson equations that consume sparse grouped `Z` directly.
Variance estimation still uses the dense covariance-component reference path.
Each term currently has one scalar variance; intercept and slope
terms in this facade are therefore independent. Use
`CorrelatedLinearMixedModel` when within-group covariance must be estimated.

`SparseLinearMixedModel` is the scalable identity-residual alternative for
independent terms. It evaluates ML/REML through sparse random-coefficient
precision equations, uses reusable minimum-degree sparse Cholesky, and avoids
forming observation-scale covariance matrices. The dense path remains the
reference for Satterthwaite/Kenward-Roger and full PEV calculations.

`CorrelatedLinearMixedModel` adds unstructured grouped intercept/slope blocks
using a positive-definite Cholesky covariance parameterization. Its current
likelihood evaluator is dense; moving correlated blocks onto the sparse
equation optimizer remains performance work.

## Pedigree animal-model REML

`Pedigree` accepts individuals in any order, validates parent references and
ancestry cycles, and constructs the numerator relationship matrix with
inbreeding by the tabular method. A `null` parent is unknown and is assumed to
be an unrelated, noninbred founder source.

```java
Pedigree pedigree = Pedigree.of(List.of(
    PedigreeIndividual.founder("sire"),
    PedigreeIndividual.founder("dam"),
    new PedigreeIndividual("offspring", "sire", "dam")
));

PedigreeRemlResult fit = PedigreeReml.fit(
    y, x, observationAnimalIds, pedigree);
```

The fitted animal model is

```text
y = X beta + Z a + e
a ~ N(0, additiveVariance A)
e ~ N(0, residualVariance I)
```

where observation identifiers define `Z`. They may repeat, and ancestors with
no phenotype may remain in the pedigree. The result includes additive and
residual variances, observation-scale narrow-sense heritability, BLUP breeding
values for every pedigree individual, prediction-error variances, and
reliabilities. Variance estimation currently forms dense `A` and `Z A Z'`;
afterward, BLUPs, PEVs, and reliabilities are obtained from Henderson equations
using sparse pedigree `A^-1` directly, including unphenotyped ancestors.

`Pedigree.sparseRelationshipMatrixInverse()` provides the `pedigreemm::getAInv`
building block directly in compressed-sparse-row form. It uses the pedigree
`T D T'` structure and does not invert dense `A`.

`SparsePedigreeReml` also uses that `A^-1` directly during variance estimation,
so neither dense `A` nor `Z A Z'` is required. The original `PedigreeReml`
remains the dense reference and currently provides richer PEV/reliability
output.

## Time-series models and LMMs with ARIMA errors

`Arima` provides Gaussian AR, AR(n), MA, ARMA, ARIMA, and multiplicative
seasonal ARIMA models. Coefficient signs match R: positive MA coefficients enter
the observation equation with a positive sign. AR and MA parameters are
transformed during optimization to preserve stationarity and invertibility.

```java
ArimaResult ar = Arima.fit(y, ArimaOrder.ar(1));
ArimaResult arma = Arima.fit(y, ArimaOrder.arma(2, 1));

ArimaOptions seasonal = ArimaOptions.builder()
    .seasonalOrder(SeasonalArimaOrder.of(1, 1, 1, 12))
    .build();
ArimaResult sarima = Arima.fit(
    y, ArimaOrder.arima(1, 1, 1), seasonal);

ArimaForecast nextYear = sarima.forecast(12, 0.95);
```

Results include the ordinary and seasonal coefficients, mean or drift,
innovations and their variance, conditional Gaussian log likelihood, AIC,
AICc, BIC, convergence metadata, and forecasts with innovation-based normal
intervals. Residual ACF and Ljung-Box diagnostics are available directly from
the result; `TimeSeriesDiagnostics` also exposes ACF and Durbin-Levinson PACF.
`includeMean` applies to undifferenced models. `includeDrift` is an explicit
option for models with exactly one ordinary or seasonal difference.

The first fitter uses conditional sum of squares/conditional Gaussian
likelihood. `ExactArma` is a separate exact stationary Gaussian estimator based
on the full Toeplitz covariance. It supports missing observations by exact
Gaussian marginalization, coefficient covariance/SE, and shared parameters
across independent panel series. `AutomaticArima.select` performs exhaustive
small-order AICc selection. Diffuse state-space likelihood for integrated or
seasonal models remains distinct from the exact stationary estimator.

`ArimaErrorLinearMixedModel` profiles the ordinary mixed-model variance
components by REML at each stationary ARMA/SARMA error-correlation candidate:

```java
ArimaErrorLmmResult fit = ArimaErrorLinearMixedModel.fit(
    y, fixedDesign, randomTerms, ArimaOrder.arma(1, 1),
    ArimaErrorLmmOptions.builder()
        .remlOptions(RemlOptions.builder()
            .degreesOfFreedomMethod(
                DegreesOfFreedomMethod.KENWARD_ROGER)
            .build())
        .build(),
    BackendPolicy.PREFERRED);
```

For ARIMA errors with differencing, the response, every fixed-effect column,
and every random-effect design column are differenced together before the
stationary ARMA-error REML fit. The returned LMM is therefore on the differenced
scale. A level intercept becomes zero under differencing and is rejected;
represent drift with a level-scale time-trend column. The ARIMA-error LMM API
still treats all rows as one ordered series; `ExactArma.fitPanel` supplies
block-independent shared ARMA estimation outside the LMM.

## Generalized mixed models: PQL with REML working fits

For binomial and Poisson mixed models, `GlmmPql` performs first-order penalized
quasi-likelihood iterations. Each iteration constructs the GLM working response
and its known diagonal residual covariance, fits the random covariance scales
by REML, and updates the observation-scale random-effect BLUP:

```java
GlmmPqlResult fit = GlmmPql.fit(
    y, x, GlmFamilies.poisson(),
    List.of(new VarianceComponent("group", observations, groupRelationship))
);
```

The random covariance bases use the same observation-scale `K` representation
as dense REML. `GlmmPqlResult` reports fixed effects, random covariance scales,
the random and total linear predictors, conditional fitted means/deviance,
outer convergence, the final working REML likelihood, and backend provenance.

PQL is an approximation, not an exact non-Gaussian restricted likelihood. It
can be biased for binary data with small clusters, rare outcomes, or large
random-effect variance. JLinAlg deliberately rejects Gaussian input in this API;
Gaussian mixed models should use exact `Reml` directly. Laplace or adaptive
quadrature GLMM likelihoods are future, separately named estimators.

The same PQL estimator is available with a numerator-relationship random
effect through `PedigreeGlmmPql`:

```java
PedigreeGlmmPqlResult fit = PedigreeGlmmPql.fit(
    y, x, GlmFamilies.binomial(), observationAnimalIds, pedigree);
```

REML, pedigree REML, GLMM PQL, and pedigree GLMM PQL all return beta, SE,
t statistics, DF, and p-values through `associationStatistics()`. The PQL
values are working-linear-model approximations rather than exact marginal
GLMM likelihood inference. Satterthwaite and Kenward-Roger options propagate to
the final PQL working REML fit; for non-Gaussian models both remain PQL-based
approximations.

## lme4 and pedigreemm compatibility direction

The existing dense APIs are compatibility-stable reference paths. Current
convenience support covers independent random intercepts/slopes, arbitrary
`Z Z'` covariance bases, multiple crossed/nested terms, conditional modes/PEVs,
pedigree Gaussian REML, pedigree binomial/Poisson/negative-binomial PQL,
profile ML, compiled fixed/mixed formulas, correlated random blocks, sparse
independent-term variance estimation, ML model comparison, and singular-fit
diagnostics. Marginal/conditional response simulation and deterministic,
parallel parametric bootstrap inference are available for ordinary Gaussian
LMMs and pedigree animal models. It does not yet implement the full `lme4`
surface: sparse correlated-block estimation, Laplace/AGQ GLMMs, and profile-
likelihood intervals remain staged work. `PreparedLinearMixedModel` retains
a model structure for warm-started response refits. `MixedModelPrediction`
provides marginal or conditional new-data prediction with an explicit
`allowNewLevels` policy. `SparsePedigreeMixedModel` composes multiple pedigree
precision terms and ordinary independent random terms in one sparse fit.

Grouped `Z` and pedigree `A^-1` now feed sparse precision equations so cost
tracks random coefficients and nonzeros rather than an observation-scale dense
`n x n` covariance. Dense execution remains valuable for numerical checks,
finite-DF adjustments, correlated blocks, and workloads where JDistlib GPU
routing is advantageous.

## Compute policy

The default `BackendPolicy.PREFERRED` order is:

1. GPU through JDistlib automatic workload routing;
2. oneMKL;
3. OpenBLAS;
4. the portable Java CPU backend.

Strict CUDA, OpenCL, Vulkan, oneMKL, OpenBLAS, GPU, automatic, and CPU policies
are also available. Results retain the selected backend and device description.
FP64 is used throughout.

## Compiled formulas and high-throughput association

`Formula.compile()` supports fixed effects, treatment or sum contrasts,
interactions, offsets, and weights. `MixedFormula.compile()` additionally
supports sparse `(1|group)` and independent `(0+x|group)` terms. Compilation
builds contiguous matrices once; fitting does not parse formulas.
It also expands `(1+x||group)` into independent terms, supports nested grouping
such as `(1|site/subject)`, and compiles `(1+x|group)` to a true
Cholesky-parameterized correlated block.

`RemlAssociationScanner` implements a batched P3D/EMMAX-style GWAS/TWAS path.
It estimates a null model once, caches its GLS projection, mean-imputes missing
dosages when requested, and evaluates marker blocks using the selected
JDistlib matrix backend. It returns marker beta, SE, t statistic, and p-value
without refitting variance components per marker. Marker blocks can also be
scanned with caller-selected CPU parallelism.

`ParallelAssociationEngine` provides ordered, bounded-parallel exact refits
when either one predictor or the response changes. Built-in adapters cover
OLS, GLM, REML/LMM, pedigree REML, GLMM/PQL, pedigree GLMM/PQL, ridge, and
ARIMA-error LMM; the `AssociationFitter` interface accepts any additional model
that produces coefficient association statistics. `FastOlsAssociation` avoids
per-variable refits by factoring the shared covariates once and residualizing
predictors in accelerated blocks. Association outputs include effect size,
SE, t/z statistic, regular p-value, `log10(p)`, and GWAS-standard `-log10(p)`;
the log forms use direct distribution log tails to avoid ordinary p-value
underflow. See [the association engine guide](docs/association-engine.md)
for execution-policy guidance and examples.

`FastGlmAssociation` fits one covariate-only GLM, caches its final IRLS metric,
and performs block efficient-score tests with one-step effect estimates. It is
the fast first-stage binary/count scan; exact alternative-model refits remain
available through `ParallelAssociationEngine`.

The file-backed pipeline accepts CSV, TSV, VCF, VCF.gz, BCF, and BGEN rather
than PLINK. It aligns samples before computing MAC/MAF/missingness/quality QC,
streams genotype blocks into fast OLS or frozen-null REML, and can write tested,
excluded, and failed rows incrementally as CSV/TSV. Feature-by-sample omics
matrices support Winsorization, log, z-score, and rank inverse-normal transforms
for TWAS/EWAS/PWAS. Continuous-trait Burden, SKAT, and SKAT-O work with either
an unrelated OLS null or the same retained related-sample REML projection.
See [the pipeline guide](docs/gwas-twas-pipeline.md) for formats, examples, and
explicit current boundaries.

## Mendelian randomization

The first MR layer accepts validated biallelic SNP associations and explicitly
harmonizes outcome effects to the exposure effect allele:

```java
HarmonizationResult harmonized = AlleleHarmonizer.harmonize(
    exposureAssociations, outcomeAssociations);

MrAnalysisResult mr = MendelianRandomization.analyze(
    harmonized.instruments());
```

Harmonization handles allele swaps and strand complements. Palindromic SNPs
require effect-allele frequencies that select exactly one orientation within a
configurable tolerance; otherwise they are excluded. Missing, duplicate,
mismatched, ambiguous, frequency-inconsistent, and zero-exposure-effect
variants have structured exclusion reasons rather than disappearing silently.

The uncorrelated-instrument analysis reports first-order Wald ratios, fixed and
multiplicative-random-effect IVW, MR-Egger slope and intercept, a weighted
median with reproducible parametric-bootstrap inference, approximate per-SNP
and mean F statistics, Cochran Q / I-squared, I-squared GX, leave-one-out IVW,
and data-dependent warnings. At least three instruments are required for the
combined analysis because MR-Egger and leave-one-out are included.

For variants retained in LD, generalized IVW and generalized MR-Egger accept a
positive-definite LD correlation matrix in exactly the harmonized instrument
order:

```java
CorrelatedMrEstimate generalized = CorrelatedMendelianRandomization.ivw(
    harmonized.instruments(), alleleAlignedLd, true);
```

LD correlation signs must correspond to each instrument's reported exposure
effect allele. MR-Egger reorientation updates the signs internally. Generalized
fits record backend provenance and use the same accelerated compute policy as
OLS and REML.

Additional sensitivity modules provide effect/SE/sample-size Steiger filtering,
Huber robust adjusted-profile-score MR, analytic PRESSO-style global/outlier
diagnostics, contamination-mixture grid profiling, multivariable IVW and
multivariable Egger, overlap-aware errors-in-variables IVW, and
selection-adjusted winner's-curse correction. These are sensitivity estimators,
not automatic proof of causality. Instrument selection, phenotype-scale
assumptions, and the exclusion restriction remain scientific responsibilities
of the caller.

## Meta-analysis and meta-regression

`MetaAnalysis` pools study-level effects with positive sampling SEs. It supports
inverse-variance fixed effects and random effects with REML, DerSimonian-Laird,
or Paule-Mandel estimates of tau-squared. Results include the pooled effect,
SE, z/t statistic, regular/log10 p-values, confidence and random-effects
prediction intervals, Cochran Q and its p-value, tau/tau-squared, I-squared,
H-squared, normalized study weights, and backend provenance.

```java
MetaAnalysisResult pooled = MetaAnalysis.fit(studies,
    MetaAnalysisOptions.builder()
        .tauSquaredEstimator(TauSquaredEstimator.REML)
        .inferenceMethod(MetaInferenceMethod.HARTUNG_KNAPP)
        .build(),
    BackendPolicy.PREFERRED);
```

`MetaRegression` accepts one or more numeric moderator columns, optionally adds
an intercept, and uses the same fixed/random heterogeneity estimators. It
returns named coefficient effect sizes, SE/statistic/p-values (including log10
forms), the full coefficient covariance, residual `Q_E`, omnibus moderator
`Q_M`, residual I-squared/H-squared, heterogeneity R-squared, and weights.

## SuSiE fine mapping

`Susie` implements iterative Bayesian stepwise selection for individual data,
standardized z-score/LD summary data, or caller-supplied `X'X`, `X'y`, and
`y'y` sufficient statistics. It reports posterior means, PIPs, per-effect
posterior inclusion matrices, residual variance, convergence, and credible
sets with LD purity. Individual inputs are centered/scaled once and returned
on their original coefficient scale. Summary LD is checked for correlation
structure and positive semidefiniteness before iteration.

## Structural equation models

`SemModel` specifies observed-variable directed paths, variances, covariances,
fixed values, and equality constraints through shared labels. `Sem` fits the
model by Gaussian covariance-structure maximum likelihood and reports parameter
SE/z/p, implied covariance, likelihood chi-square, CFI, TLI, RMSEA, SRMR,
AIC, and BIC. Variances use a positive log parameterization. Complete-case
omission is optional. Latent-variable measurement models, ordinal thresholds,
robust sandwich corrections, and full-information missing-data likelihood are
future extensions rather than being silently approximated.

## Numerical scope

See [docs/numerical-contract.md](docs/numerical-contract.md) for assumptions,
likelihood definitions, convergence behavior, and current limitations.

Cross-language regression fixtures compare JLinAlg against base R, `nlme`,
`rrBLUP`, and `MendelianRandomization`, including base-R inverse-variance
meta-analysis calculations. See
[docs/r-reference-validation.md](docs/r-reference-validation.md) for exact
versions, cases, tolerances, and the reproducible R generator.

## License

JLinAlg is licensed under the GNU General Public License, version 2 or later
(`GPL-2.0-or-later`).
