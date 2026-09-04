# MR, time-series, SuSiE, and SEM scope

> **v0.1.0 performance status:** MR and SEM remain initial performance paths.
> SuSiE is directly validated against susieR and benchmarked on the package's
> official `N3finemapping` vignette data.

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

## Time series

`Arima` remains the fast conditional estimator and supports integration,
seasonality, forecasting, and ARIMA-error LMMs. `ExactArma` evaluates the full
stationary Gaussian Toeplitz covariance. Missing values are marginalized by
subsetting that covariance, and independent panel likelihoods are summed under
shared coefficients. This is exact for stationary ARMA, but is not an exact
diffuse Kalman likelihood for integrated models.

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

`SemModel` uses the observed-variable RAM relationship
`Sigma = (I - A)^-1 S (I - A)^-T`. Free residual variances are optimized on
the log scale. Shared labels impose equality constraints. Current estimation
uses complete-case covariance ML; latent measurement variables, mean
structures, ordinal likelihoods, robust corrections, modification indices,
and FIML missingness are outside this first engine.
