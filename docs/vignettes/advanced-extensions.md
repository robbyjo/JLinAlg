# Advanced extensions: methods and validation

The September 8 audit replaced several incomplete implementations, not just
their descriptions. The [accuracy and performance report](../advanced-validation.md)
links independent R fixtures, measured timings, and reproduction commands.
No timings imply complete R-package parity or a universal speed advantage.

## SEM

`SemModel` supports latent measurement paths and structural intercepts.
`Sem.fit` jointly estimates the RAM model; `SemFiml.fit` maximizes that same
model's observed-pattern Gaussian likelihood directly. It does not fit a
saturated covariance and treat a subsequent complete-data fit as FIML.
`SemOrdinal.fit` jointly estimates ordinal thresholds and structural parameters
using a probit pairwise likelihood. Full parameter covariance, case/cluster
sandwich corrections, continuous-model efficient modification indices, and
indirect-effect delta-method inference are available.

Checked-in `lavaan` references cover latent paths, means, FIML, robust ML, and
ordinal PML. Ordinal PML is not DWLS/WLSMV. Missing/mixed ordinal models,
multigroup SEM, ordinal modification indices, and FIML robust scaled fit
statistics remain unsupported. See [SEM](sem.md) for identification and examples.

## Meta-analysis

`MetaMultilevelRegression` jointly estimates ML/REML random moderator covariance
blocks over a full sampling covariance matrix. `MetaClusterRobust` supports
CR0, CR1, CR1P, CR1S, and CR2, including coefficient-specific Satterthwaite
inference for CR2. `MetaBiasCorrections` provides PET/PEESE and iterative
Duval–Tweedie L0/R0/Q0 trim-and-fill, replacing the former median reflection.
Effect-size construction includes stable large-sample Hedges corrections.

`metafor` and `clubSandwich` fixtures gate covariance, likelihood, bias
corrections, and inference. Computation is dense; covariance optimization is
local, and variance-component profile intervals are not provided. Bias
diagnostics do not establish or remove publication bias. See [meta-analysis](meta-analysis.md).

## Mixed models

The sparse path estimates unstructured random-effect covariance with prepared
sparse equation refits. Formula compilation supplies the design, weights,
offsets, and grouping terms; profile refits hold the target fixed while
reoptimizing nuisance parameters. Finite-DF inference uses sparse likelihood
and fixed-covariance derivatives, not an observation-sized covariance inverse.
The reference suite compares likelihood, coefficients, covariance, finite DF,
and profile endpoints to `lme4`, `lmerTest`, and `pbkrtest`.

See [the mixed-model roadmap](../lme4-pedigreemm-roadmap.md) and
[linear models](linear-models-and-glms.md). Boundary or singular covariance
fits do not justify ordinary interior Wald/finite-DF approximations.

## Time series

`DiffuseArima` uses exact diffuse Kalman recursions for integrated and seasonal
models, including missing observations. `SparseMissingSeries` uses
state-sized filtering rather than an observation-by-observation dense
covariance. Forecasts propagate the filtered state, not a zero-filled series.
References compare to `stats::arima`; a bounded-heap test covers 50,000 positions.

R's finite diffuse initialization has a small, documented likelihood difference
from symbolic exact diffuse initialization. Diffuse coefficient covariance,
general regression terms, and historical smoothing are not exposed.
See [time series](time-series.md).

## MR workflow

The `mr-estimate` CLI exposes generalized LD-aware IVW/Egger, overlap-aware
IVW, native SVG plots, and actual conditional association calculations.
`conditional` requires an explicit `--joint` or `--condition-on` choice;
`secondary-signals` performs conditional forward/backward signal selection.
Conditional instrument p-values and causal-estimate p-values are different
quantities with separate output names.

```powershell
java -jar build/cli/jlinalg-0.3.2.jar mr-estimate `
  --input harmonized.tsv --method conditional --joint --ld ld-correlation.tsv
```

Signed, allele-aligned LD and the documented Gaussian summary-score assumptions
are required. Generalized MR is not itself a conditional GWAS test, and the
older marginal LD clumper remains explicitly marginal.
See [end-to-end MR](mr-end-to-end.md).

## GLMM

`GlmmQuadrature` adapts Gauss–Hermite nodes to each group's random-intercept
mode and curvature, refines the node count, checks the zero-variance fit,
and exposes numerical diagnostics and joint covariance. Binomial trials,
Poisson responses, and offsets are supported. Independent `integrate` and
`lme4` high-node AGQ fixtures gate likelihood and inference.

This is a scalar independent Gaussian random-intercept alternative, not
multidimensional or pedigree AGQ. Refinement differences are diagnostics, not
certified integration error bounds. See [GLMM and pedigree](pedigree-and-glmm.md).

## Penalized and other regression families

`SelectionAwarePenalizedInference` preserves weights and penalty factors while
selecting on one sample and performing inference on independent held-out rows.
`PolyhedralSelectiveInference` additionally conditions on the active set and
signs at a fixed LASSO/elastic-net penalty, using truncated-normal p-values and
inverted confidence intervals. Its Gaussian noise SD must be known independently;
response-tuned lambda and dependent selection/inference rows invalidate these
contracts. `glmnet`/`selectiveInference` fixtures check the LASSO case, with
analytic elastic-net and null-calibration checks.

Multivariate OLS reuses one QR across outcomes. Multinomial and smoothed-quantile
fits use score-checked optimizers; `QuantileRegression.fitExact` adds a
nonsmoothed primal-dual LP with optimality certificates. Supersmoothing implements R's running-window
algorithm, including weights, periodicity, fixed spans, and bass control.
Kernel prediction remains finite far outside the sample, and the partially
linear API adds Robinson residualization and slope HC3 inference.
See [regression families](regression-families.md) and
[selection inference](linear-models-and-glms.md).

## Run the integration gate

```powershell
.\gradlew.bat check benchmarkClasses
```

The R generators are checked in; tests consume frozen fixtures and do not
require R or download packages at test time. Performance comparisons have
explicit workload, initialization, and inference-work caveats in the report.
