# MR, time-series, SuSiE, and SEM scope

> **Performance status:** SEM is directly validated against `lavaan` and
> benchmarked on a TOPMed cardiometabolic path model and joint latent/ordinal/FIML
> fixtures. MR conditional/generalized/overlap-aware paths now have paired R timings.
> See the [September 8 audit](advanced-validation.md). SuSiE is validated against susieR and benchmarked
> on the package's official `N3finemapping` vignette data.

## Mendelian randomization

The core independent- and LD-aware estimators are supplemented by:

- `SteigerFiltering` for instrument and aggregate directionality;
- `RobustMendelianRandomization.raps` for Huber adjusted-profile scoring;
- `MrPresso` for fast analytic, robust-center outlier diagnostics;
- `ContaminationMixture` for valid/invalid-instrument mixture profiling;
- `MultivariableMendelianRandomization` for IVW and Egger direct effects;
- `OverlapAwareMendelianRandomization` for per-instrument sampling covariance;
- `WinnerCurseCorrection` for selection-adjusted normal likelihood.

The PRESSO implementation is explicitly analytic and does not claim the
simulation calibration of the R MR-PRESSO package. Exposure/outcome covariance
must be supplied by the caller when samples overlap.

`ConditionalAssociation` fits a Gaussian summary-score model with signed,
allele-aligned LD. Its conditional instrument p-values are not causal MR
p-values. `SecondarySignalClumper.select` implements conditional forward/backward
selection; the older `clump` method is explicitly marginal LD pruning. The
`mr-estimate` CLI and native SVG output expose these separate workflows; see
[the end-to-end vignette](vignettes/mr-end-to-end.md).

## Time series

`Arima` remains the fast conditional estimator and supports integration,
seasonality, forecasting, and ARIMA-error LMMs. `ExactArma` evaluates stationary
Gaussian likelihood with state-sized missing-observation filtering; independent
panel likelihoods are summed under shared coefficients. `DiffuseArima` adds
symbolic exact diffuse recursions for integrated and seasonal models, including
missing values and filtered-state forecasts. R's finite initialization and
Java's exact diffuse normalization are distinguished in the
[time-series vignette](vignettes/time-series.md).

## SuSiE

`Susie` implements IBSS with a Gaussian single-effect prior. Summary inputs
use susieR's finite-sample transformation: for
`a_j = (n - 1)/(z_j^2 + n - 2)`, `X'X = (n - 1) R`,
`X'y_j = sqrt((n - 1) a_j) z_j`, and `y'y = n - 1`. For other scales, use the
sufficient-statistics entry point. Credible-set purity is the minimum absolute
LD inside each set. `priorVariance` is the absolute variance on the standardized
predictor scale; to reproduce it in susieR, pass
`scaled_prior_variance = priorVariance / var(y)` and disable prior-variance
estimation.

## SEM

`SemModel` uses RAM covariance `Sigma = (I-A)^-1 S (I-A)^-T` and means
`mu = (I-A)^-1 intercept`, selecting observed margins from the joint latent
distribution. Log variances and shared-label equality constraints are supported.
`Sem` jointly fits latent/mean models; `SemFiml` directly maximizes the
constrained observed-pattern likelihood. `SemOrdinal` jointly fits thresholds
and RAM parameters with ordinal probit pairwise likelihood. Robust/cluster
covariance, continuous-model efficient modification indices, and indirect delta
inference are explicit APIs. Ordinal PML is not DWLS/WLSMV. See the
[SEM vignette](vignettes/sem.md) for identification and unsupported extensions,
and [TOPMed performance](topmed-sem-performance.md) for the earlier observed model.
