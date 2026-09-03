# MR, time-series, SuSiE, and SEM scope

> **v0.1.0 performance status:** MR, SuSiE, and SEM are implemented and tested,
> but have not yet been optimized or benchmarked for large or high-throughput
> workloads. Their APIs should be treated as initial in this release.

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
assume standardized quantitative traits and predictors: `X'X = n R`,
`X'y = sqrt(n) z`, and `y'y = n`. For other scales, callers should use the
sufficient-statistics entry point. Credible-set purity is the minimum absolute
LD inside each set.

## SEM

`SemModel` uses the observed-variable RAM relationship
`Sigma = (I - A)^-1 S (I - A)^-T`. Free residual variances are optimized on
the log scale. Shared labels impose equality constraints. Current estimation
uses complete-case covariance ML; latent measurement variables, mean
structures, ordinal likelihoods, robust corrections, modification indices,
and FIML missingness are outside this first engine.
