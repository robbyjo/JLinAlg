# Advanced extensions

This vignette records the bounded implementations that close the remaining
analysis-extension backlog. They are deterministic APIs with explicit
assumptions; they do not silently claim full package-level parity with R.

## SEM

The SEM package now includes:

- `LatentMeasurement`, a deterministic principal-factor measurement model for
  continuous indicators;
- `SemMeanStructure`, which reports intercept/mean estimates and their
  covariance;
- `SemInference.ordinalThresholds`, which estimates normal-response thresholds
  from ordinal category counts;
- `SemFiml`, an EM saturated Gaussian mean/covariance estimator that preserves
  missing-data patterns before passing the covariance to the structural fit;
- `SemInference.sandwich` for case/cluster score rows and a caller-supplied
  inverse bread;
- `SemInference.modificationIndices` for score/information diagnostics; and
- `SemInference.indirect` for product-of-coefficients delta-method intervals.

The latent path is a principal-factor measurement fit rather than a complete
latent-variable RAM optimizer, and ordinal thresholds are marginal normal
thresholds rather than a polychoric ordinal likelihood. These boundaries are
returned in the API documentation so results are not mistaken for lavaan
DWLS/FIML parity.

## Meta-analysis

`MetaMultilevelRegression` fits a moderator design against a full study-by-study
sampling covariance matrix with GLS. This supports correlated effects and
multilevel covariance blocks without flattening them into independent studies.
`MetaBiasCorrections` adds deterministic PET, PEESE, and a transparent
median-reflection trim-and-fill diagnostic. The corrected effect and number of
imputed studies are returned; the trim-and-fill method is labelled explicitly
because it is not a replacement for a prespecified publication model.

## Mixed models

`SparseUnstructuredCorrelatedModel` estimates a lower-Cholesky random-effect
covariance shape through deterministic sparse equation refits. A single
correlated block can also be reached from `CompiledMixedFormula` through
`fitSparseUnstructured`, preserving formula compilation and sparse factorization
without materializing an observation-scale covariance. Existing
`MixedModelProfile` supplies profile-likelihood crossings for formula or model
refit callbacks.

## Time series

`DiffuseArima` exposes integrated ARIMA fitting with explicit diffuse-state
counting and effective likelihood observations. `SparseMissingSeries` routes
stationary ARMA fits with arbitrary missing positions through the exact observed
pattern covariance path in `ExactArma`; contiguous complete series retain the
Toeplitz innovation path.

## MR workflow

`mr-estimate` now accepts `--ld MATRIX` for generalized fixed/random IVW,
generalized Egger, and the `conditional` alias, and `--sampling-covariance FILE`
for overlap-aware IVW. Generalized estimators compute conditional p-values
against the supplied LD covariance rather than treating correlated instruments
as independent. Native SVG output remains available through `--plot`.

Example:

```powershell
java -jar build/cli/jlinalg-0.2.0.jar mr-estimate `
  --input harmonized.tsv --method ivw-generalized-random `
  --ld ld-correlation.tsv --plot mr.svg
```

## GLMM

`GlmmQuadrature` provides a fixed 10-node Gauss--Hermite marginal likelihood
for binomial or Poisson-style one-dimensional random-intercept models. It is a
deterministic quadrature alternative to PQL and first-order Laplace when the
random effect is scalar. Higher-dimensional random effects still require the
sparse Laplace/PQL paths or an external adaptive-quadrature reference.

## Penalized regression

`SelectionAwarePenalizedInference` implements deterministic data-splitting:
the first fraction of rows selects the active set with LASSO/elastic net, and
the held-out rows provide the OLS inference fit. This produces inference that
is independent of the selection sample, unlike ordinary post-selection OLS;
the result records both sample sizes and the selected indices.

## Verification

Focused coverage is in `SemExtensionsTest`, `MetaExtensionsTest`,
`MixedExtensionsTest`, `TimeSeriesExtensionsTest`, `MrEstimatorExtensionTest`,
`GlmmQuadratureTest`, and `SelectionAwareInferenceTest`. Run the complete
repository check before publishing:

```powershell
.\gradlew.bat check
```
