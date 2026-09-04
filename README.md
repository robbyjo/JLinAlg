# JLinAlg

**[Explore the project website](https://robbyjo.github.io/JLinAlg/)** for the
complete feature map and browser-friendly worked vignettes.

JLinAlg implements Java linear and mixed-model algorithms on top of
[JDistlib 0.10.1](https://github.com/robbyjo/JDistlib/releases/tag/v0.10.1).
Version 0.1.0 provides ordinary least squares (OLS), generalized
linear models (GLMs), dense Gaussian restricted maximum likelihood (REML),
pedigree animal-model REML, and penalized-quasi-likelihood generalized linear
mixed models (GLMM PQL). Gaussian ridge, LASSO, and elastic-net regression and
a first summary-statistics Mendelian-randomization (MR) layer are also included.

> **v0.1.0 performance status:** MR and SEM are implemented and tested, but have
> not yet received the large-workload optimization and benchmarking
> applied to JLinAlg's performance-tuned model and association paths. Treat
> their APIs as initial and profile them on representative data before
> high-throughput use.

## Requirements and build

- A JDK 17 or newer.
- Network access on the first build. The build downloads the pinned
  `jdistlib-all-0.10.1.jar` and verifies its SHA-256 digest before compiling.

On Windows:

```powershell
.\gradlew.bat check
```

On Linux or macOS:

```shell
./gradlew check
```

No system Gradle installation is needed.

## Executable command line

Build the self-contained command-line JAR with:

```powershell
.\\gradlew.bat executableJar
```

The artifact is written to `build/cli/jlinalg-<version>.jar` and includes its
runtime dependencies. A fixed-effect omics scan can then be run as:

```powershell
java -jar build/cli/jlinalg-0.1.0.jar \`
  --omics methylation.tsv \`
  --pheno phenotype.tsv \`
  --id IID \`
  --formula "trait ~ age + sex + <omics>" \`
  --transform "<omics>=mvalue(epsilon=1e-6)|zscore()" \`
  --out ewas-results.tsv
```

`<omics>` is the one formula term sourced from the omics input. Every other
response, fixed-effect, interaction, and random-effect name must be a phenotype
column. Without `--omics`, the command fits one phenotype-only model and emits
all fixed-effect coefficients. Gaussian fixed-only formulas resolve to OLS;
Gaussian formulas containing terms such as `(1|Batch)` resolve to REML. A
non-Gaussian `--family` resolves to GLM or PQL GLMM, and `Surv(time,event)`
uses Cox regression.

Add a genomic relationship matrix with `--grm FILE`. Its presence makes an
otherwise fixed Gaussian or non-Gaussian formula resolve to LMM or GLMM,
respectively; a Cox formula uses the GRM as a Gaussian kinship frailty. For
repeated observations, `--individual-id COLUMN` names the phenotype column
that maps rows to GRM individuals. It defaults to the `--id` column.

```powershell
java -jar build/cli/jlinalg-0.1.0.jar `
  --pheno phenotype.tsv --id observation_id --individual-id IID `
  --formula "trait ~ age + sex" --grm cohort `
  --out trait-grm.tsv
```

`--grm` accepts either a labeled square CSV/TSV (first column contains row
IDs; remaining headers contain column IDs) or a GCTA prefix backed by
`PREFIX.grm.bin` and `PREFIX.grm.id`. GCTA binary values are read as the
standard little-endian, single-precision lower triangle. If IIDs are not
unique, the reader exposes them as `FID:IID`. The GRM is necessarily retained
as a dense sample-by-sample covariance matrix, while the much larger omics
matrix remains block streamed. GRM provenance, format, matching column, and
memory footprint are recorded in the run log and manifest.

Delimited omics inputs have features in rows and sample IDs in the header.
VCF, BCF, and BGEN inputs use the existing streaming genotype readers.
CSV/TSV feature identifiers are inspected conservatively to infer GWAS, EWAS,
or Ensembl expression schemas; `--omics-type` overrides the inference.

The omics matrix is never materialized in full. Readers validate in constant
memory and process an adaptive block sized from current JVM heap headroom.
`--block-size N` overrides the automatic choice. Results are spooled as they
arrive, and Benjamini-Hochberg adjustment uses bounded external sorting rather
than retaining every result or p-value in memory.

The common result table contains beta, standard error, test statistic,
denominator DF, partial R-squared where defined, p-value, and BH FDR. Genotype
results additionally report ALT-oriented effects, allele frequencies, MAC,
missingness, and all-sample HWE. Recognized binary responses add separate case
and control HWE columns. Non-WGS annotations can be joined with:

```text
--annot annotation.tsv --annot-id probe_id --annot-cols chr,start,gene,strand
```

Logging is on by default. `--out results.tsv` creates `results.tsv.log` and
`results.tsv.manifest.json`; `--log` overrides the human-readable log path.
Use `--dry-run` to validate and resolve an analysis without fitting it.

Built-in transform pipelines include winsorization, shifted/log1p transforms,
z-scoring, Blom inverse-normal ranks, and EWAS M values. Trusted custom Java
transforms implement `OmicsTransformProvider`, register it through Java's
service-provider mechanism, and are loaded with `--transform-plugin FILE.jar`.

## Vignettes and worked examples

The [vignette index](docs/vignettes/README.md) provides end-to-end examples for
every implemented feature group: OLS/GLM and penalized regression, REML/LMM,
GAM/GAMM and distributional/vector additive models, pedigree and GLMM
PQL/Laplace, formulas/backends, association and GWAS/TWAS, Mendelian
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
Positive observation weights and nonnegative per-coefficient penalty factors
are supported; a zero factor leaves a covariate unpenalized.

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

PenalizedRegression.Prepared prepared =
    PenalizedRegression.prepare(y, x, options);
PenalizedRegressionPath ridgePath =
    prepared.automaticPath(100, 1e-4, 0.0);
PenalizedRegressionPath lassoPath =
    prepared.automaticPath(100, 1e-4, 1.0);
```

`automaticPath` creates a descending log-spaced lambda path and uses warm
starts. `PenalizedRegressionCrossValidation` performs reproducible K-fold
selection and returns both the minimum-error lambda and the conservative
one-standard-error choice. Preprocessing is learned separately inside each
training fold. A prepared design reuses centering, scaling, response products,
and the covariance matrix across alpha values or user-supplied paths.

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
Gaussian mixed models should use exact `Reml` directly.

`SparseGlmmLaplace` is the marginal-likelihood alternative for grouped and
pedigree models. It consumes sparse `RandomEffectTerm` designs and
`SparsePrecisionMatrix` coefficient precisions directly. Its `Prepared` API
owns the selected backend for an association scan and reuses one symbolic and
numeric sparse Cholesky factor per worker, avoiding both observation-scale
covariance matrices and repeated backend discovery. `GlmmLaplace` remains the
dense reference implementation; adaptive quadrature is not implemented.

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

## Cox proportional-hazards and frailty models

`CoxRegression` fits right-censored or counting-process Cox models with
strata, offsets, Breslow or Efron ties, hazard ratios, SE/z/p inference, and a
stratum-specific baseline cumulative hazard. The common right-censored path
sorts each stratum once and accumulates risk-set moments in descending time;
delayed-entry data use the general start-stop reference path.

```java
CoxSurvivalData survival = CoxSurvivalData.rightCensored(time, event);
CoxResult fit = CoxRegression.fit(survival, covariates);

double hazardRatio = fit.hazardRatios()[0];
double se = fit.standardErrors()[0];
double p = fit.pValues()[0];
```

`CoxMixedModel` adds one or more Gaussian log-frailty terms and profiles their
variances with a Laplace-approximated marginal partial likelihood. Independent
grouped terms can be converted from `RandomEffectTerm`; arbitrary precision
matrices are accepted through `CoxRandomEffectTerm`. `CoxPedigreeFrailty`
uses the pedigree's sparse-built `A^-1` directly as the additive frailty
precision and returns named individual modes.

```java
CoxRandomEffectTerm center = CoxRandomEffectTerm.independent(
    RandomEffectTerm.randomIntercept("center", centerIds));
CoxMixedResult mixed = CoxMixedModel.fit(
    survival, covariates, List.of(center));

CoxPedigreeResult pedigreeFit = CoxPedigreeFrailty.fit(
    survival, covariates, animalIds, pedigree);
```

The fixed model is an ordinary Cox partial-likelihood estimator. Mixed and
pedigree models are explicitly Laplace-approximated Gaussian frailty models,
not Gaussian REML and not gamma frailty. See the
[survival vignette](docs/vignettes/cox-survival.md).

## Cryptic relatedness and genomic relationship matrices

`GenomicRelationshipMatrix` constructs an additive GRM from aligned,
variant-by-sample dosages using mean imputation and per-variant
`sqrt(2 p (1-p))` standardization. MAF and call-rate filters are explicit, and
the final `Z Z'` product routes through JDistlib so GPU, oneMKL, or OpenBLAS can
handle the matrix-heavy step. `relatedPairs()` reports positive off-diagonal
relationships and the corresponding kinship coefficients.

```java
GenomicRelationshipMatrix grm =
    GenomicRelationshipMatrix.fromVariants(
        ldPrunedVariants, sampleIds,
        GenomicRelationshipOptions.defaults(),
        BackendPolicy.PREFERRED);

VarianceComponent cryptic = grm.varianceComponent("grm");
List<RelatednessPair> pairs = grm.relatedPairsByKinship(0.0884);
```

The variance component plugs directly into Gaussian REML, GLMM PQL,
`RemlAssociationScanner`, and therefore related-sample Burden, SKAT, and
SKAT-O. An overload expands `Z K Z'` for reordered or repeated observation
IDs. `CoxKinshipFrailty` uses the same GRM as a Laplace Gaussian frailty and
handles near-singular empirical matrices with an explicit relative diagonal
regularization.

```java
CoxMixedResult survivalFit = CoxKinshipFrailty.fit(
    survival, covariates, observationSampleIds, grm);
```

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

1. SuiteSparse CHOLMOD for sparse factorization, paired with the best
   available native dense CPU backend;
2. GPU through JDistlib automatic workload routing;
3. oneMKL;
4. OpenBLAS;
5. the portable Java CPU backend.

Strict CHOLMOD, CUDA, OpenCL, Vulkan, oneMKL, OpenBLAS, GPU, automatic, and CPU
policies are also available. The CLI exposes the same choices through
`--backend`; `--backend preferred` is the default. Results retain the
selected backend and device description. FP64 is used throughout.

The CHOLMOD JNI bridge pins SuiteSparse 7.12.2 by source URL and SHA-256.
It retains one symbolic analysis per worker and repeats only numeric
factorization while the sparsity pattern is unchanged. Build and stage the
host library with:

```powershell
.\gradlew.bat nativeCholmodJar
```

The build uses CMake and the system BLAS/LAPACK. On this Windows development
host it detects oneMKL; another installation can be selected with
`-PcholmodMklRoot=PATH`. An existing SuiteSparse checkout can be supplied
with `-PsuiteSparseSourceDir=PATH`. During development the loader finds
`build/native/cholmod/bin/Release/jlinalg_cholmod.dll`; packaged
platform-specific resources and an explicit
`-Djlinalg.cholmod.library=PATH` are also supported. If the native bridge is
absent, `PREFERRED` follows the fallback order above, while strict
`CHOLMOD` fails with the native loader diagnostic.

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

The v0.1.0 MR implementation is tested for numerical behavior but is not yet
optimized or benchmarked as a high-throughput pipeline.
The [end-to-end MR vignette](docs/vignettes/mr-end-to-end.md) connects public
or custom instrument preparation, LD clumping, harmonization, estimation,
diagnostics, plotting exports, bidirectional MR, and molecular-trait designs.
The [parallel xWAS MR vignette](docs/vignettes/xwas-mr-pipeline.md) adds a
bounded two-stage exposure-by-phenotype scan that reuses clumped instruments
and reserves full MR diagnostics for threshold-passing pairs.
The file-to-file CLI and reproducible Java/R results are documented in the
[xWAS MR performance report](docs/xwas-mr-cli-performance.md).

Candidate instruments can be discovered by trait in the public NHGRI-EBI GWAS
Catalog and streamed directly into MRInstruments-compatible columns. Common
trait acronyms such as `BMI` are expanded before the ontology-backed search:

    java -jar jlinalg-<version>.jar mr-instruments search \
      --trait BMI --limit 20
    java -jar jlinalg-<version>.jar mr-instruments download \
      --study GCST... --out bmi-instruments.tsv --p-threshold 5e-8

The download command reads the Catalog's standardized summary-statistics file
as a stream, keeps candidate rows at the requested p-value threshold, converts
odds or hazard ratios to log effects when necessary, and reports invalid or
non-biallelic rows. It does not claim that significance alone establishes a
valid instrument: ancestry-matched LD clumping and the usual relevance,
independence, and exclusion-restriction checks remain required.

User-provided GWAS or QTL tables can be normalized in the same way. CSV, TSV,
and gzip-compressed files are supported; common names are detected
automatically, while `--map` resolves study-specific headers:

    java -jar jlinalg-<version>.jar mr-instruments format \
      --input exposure.csv.gz --out exposure.mr.tsv --trait BMI \
      --map SNP=rsid,beta=estimate,se=stderr,eaf=frequency,\
        effect_allele=allele1,other_allele=allele2,pval=p_value

The canonical output is `Phenotype`, `SNP`, `beta`, `se`, `eaf`,
`effect_allele`, `other_allele`, `pval`, `units`, `ncase`, `ncontrol`,
`samplesize`, and `gene`. Use `--effect-scale odds-ratio` or
`--effect-scale hazard-ratio` when an explicitly mapped effect is a ratio.

Freely available LD reference databases can be listed and installed with the
CLI. Every source is normalized to the versioned, variant-major PLINK layout
described in the [LD reference format](docs/ld-reference-format.md):

    java -jar jlinalg-<version>.jar ld-db list
    java -jar jlinalg-<version>.jar ld-db download \
      --database 1000g-phase3 --location /data/ld/1000g-phase3

Clump the candidate table directly against an ancestry-matched installed panel:

    java -jar jlinalg-<version>.jar clump \
      --database /data/ld/1000g-phase3 --population EUR \
      --instrument bmi-instruments.tsv --ld-threshold 0.001 \
      --output bmi-instruments.clumped.tsv

The TwoSampleMR-compatible defaults are a 10,000 kb window, r-squared 0.001,
and index p-value threshold 1. The output preserves all input columns.

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

High-throughput scans can validate and copy row-major inputs once, then fit
all independent rows without constructing per-study objects or routing tiny
one-column systems through a general matrix backend:

```java
PreparedMetaAnalysisBatch scan = MetaAnalysis.prepareBatch(
    rowMajorEffects, rowMajorStandardErrors, features, studiesPerFeature);
MetaAnalysisBatchResult results = scan.fit(
    MetaAnalysisOptions.randomEffects(), 8);
```

The batch result is columnar and includes pooled effects, errors, test
statistics and p-values, intervals, Cochran Q, tau-squared, I-squared, and
H-squared. See [the four-cohort TOPMed profile](docs/topmed-meta-analysis-performance.md)
for the Java/R benchmark and reproduction commands.

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
structure before iteration; as in susieR's default `check_input = FALSE`, an
eager cubic positive-semidefinite factorization is not performed.

The implementation matches fixed-prior susieR results on its official
`N3finemapping` vignette data within `2e-10` and is benchmarked end to end.
On the documented i9-9900K run it took 0.0704 seconds, 4.8 times faster than
susieR 0.14.2 on the same data and host.

## SuSiE colocalization

`ColocSusie` implements the multi-signal method exposed by R
`coloc::coloc.susie`. It compares every retained pair of credible signals and
returns posterior probabilities for H0 through H4 together with the
variant-level posterior conditional on a shared causal variant. Inputs are
aligned by variant identifier, and low-posterior-overlap signal pairs are
trimmed with the same defaults as coloc.

Use `ColocSusie.analyze(firstFit, secondFit)` with JLinAlg `SusieResult`
objects. For fits produced by `susieR`, pass the selected `lbf_variable` rows
through `ColocSusieInput`; this avoids refitting and gives direct numerical
interoperability. The combination loop precomputes variant alignment and
per-signal reductions, and performs no maps or temporary allocations inside
the O(signals1 * signals2 * variants) hot path.

## Structural equation models

The v0.1.0 SEM implementation is tested but is not yet optimized or benchmarked
for large models or repeated fits.

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
