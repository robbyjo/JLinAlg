# Development TODO

Last reviewed: 2026-09-09.

This inventory tracks open work. Completed implementation and audit details are
recorded in the [advanced-method validation report](docs/advanced-validation.md),
the [v0.3.0 remaining-code audit](docs/release-0.3.0-audit.md), and the
[release notes](RELEASE_NOTES.md).

## Zelig-inspired additions — high priority

- [ ] **Unified predictions and contrasts.** Add a common API for expected
  responses, scenario differences, risk ratios, and average marginal effects
  with confidence intervals. Reuse fitted coefficient covariance and use
  analytic gradients/delta-method inference where suitable; simulation should
  be optional. Account for covariance between scenario predictions and
  distinguish population-averaged effects from predictions at average
  covariates. Keep predictive uncertainty distinct from uncertainty in the mean.
  Reference: [Zelig workflow](https://github.com/IQSS/Zelig#zelig-workflow-overview).
- [ ] **Binary probit regression.** Add a built-in binomial probit family using
  the existing `GlmFamily`/IRLS infrastructure, with stable normal-tail
  probabilities, derivatives, likelihoods, and inference.
  Reference: [Zelig probit](https://christophergandrud.github.io/Zelig/articles/zelig_probit.html).
- [ ] **Individual-level instrumental-variable regression / 2SLS.** Add general
  regression with endogenous predictors and instruments, complementing the
  existing summary-statistic Mendelian-randomization methods. Use stable QR/SVD
  projections, identification/rank checks, correct IV covariance, robust/cluster
  inference, and instrument-strength diagnostics. Do not substitute ordinary
  second-stage OLS standard errors for IV inference.
  Reference: [Zelig IV source](https://github.com/IQSS/Zelig/blob/master/R/model-ivreg.R).

## Zelig-inspired additions — medium priority

- [ ] **Rare-events logistic regression.** Implement King–Zeng coefficient-bias
  correction and population-prevalence adjustments for case-control sampling.
  Separate bias correction from prevalence correction and explicitly report
  separation or failed underlying logistic fits.
  Reference: [Zelig ReLogit source](https://github.com/IQSS/Zelig/blob/master/R/model-relogit.R).
- [ ] **Tobit / censored Gaussian regression.** Add likelihood-based regression
  for measurements censored at known limits. Expose latent-response means,
  observed-response means, and censoring probabilities separately, with stable
  tail calculations and covariance inference.
  Reference: [Zelig Tobit](https://christophergandrud.github.io/Zelig/articles/zelig_tobit.html).
- [ ] **Parametric survival regression.** Add Weibull, exponential, and lognormal
  accelerated-failure-time models with censoring, survival-time predictions,
  time ratios, and uncertainty estimates. Document distribution/scale
  parameterizations and validate against `survival::survreg`.
  References: [Zelig model catalog](https://github.com/IQSS/Zelig/tree/master/R),
  [survreg documentation](https://stat.ethz.ch/R-manual/R-devel/RHOME/library/survival/html/survreg.html).
- [ ] **Ordered logit/probit by maximum likelihood.** Add a dedicated
  cumulative-link ordinal likelihood fitter with ordered thresholds, category
  probabilities, and covariance inference. This extends the existing ordinal
  GEE, adjacent-category logits, and ordinal SEM with a distinct estimator.
  Reference: [ZeligChoice ordered logit](https://christophergandrud.github.io/Zelig/articles/zeligchoice_ologit.html).
- [ ] **Survey-design inference.** Add an explicit sampling-design abstraction
  and design-based regression inference, including sampling weights, strata,
  and primary sampling units. Ordinary regression weights alone do not supply
  survey-design variance estimates. Validate against the R `survey` package.
  Reference: [survey GLM documentation](https://github.com/cran/survey/blob/master/man/svyglm.Rd).

## Validation for new methods

Use checked-in fixtures from the underlying R estimators, independent
likelihood/derivative/covariance checks, and reproducible workload-specific
benchmarks. Zelig often wraps other packages and must not be the sole numerical
reference. Document estimator assumptions, convergence failures, and supported
scope; scenario contrasts require additional assumptions for causal interpretation.

## Remaining extensions and limits

These boundaries are explicit; this release does not claim complete parity
with every feature of lavaan, metafor, lme4, or the other R packages.

- SEM: DWLS/WLSMV, mixed/missing ordinal responses, multigroup invariance,
  ordinal modification indices, and FIML robust scaled fit statistics.
- Meta-analysis: sparse large-study covariance estimation, variance-component
  profiles, and additional structured covariance families.
- Mixed models: joint multi-DF KR tests, nonregular boundary-profile coverage,
  general correlation profiles in blocks larger than two, original-coordinate
  full PEV blocks, and joint unstructured-plus-pedigree covariance optimization.
  Formula-level pedigree mapping and automatic complete-case alignment remain open.
- Time series: diffuse coefficient covariance, general regression terms, and
  historical state smoothing.
- MR: additional generalized estimators requiring richer summary statistics;
  conditional GWAS calculations currently use the documented Gaussian score model.
- GLMM: multidimensional/crossed/pedigree adaptive quadrature beyond the new
  scalar Gaussian random-intercept path.
- Zero-inflated mixed models: a reproduced competing-mode ZINB refit is not
  certified. It now reports nonconvergence and unavailable inference instead of
  false success; see the exact case in the distributional audit. General
  multi-mode integration and global optimization remain open.
- Nonlinear models: floating-point quantization under very large response
  offsets can prevent a score certificate even when fitted SSE is near optimal.
  Such fits return nonconvergence. Mixed effects are additive Gaussian effects,
  not a general nonlinear random-parameter likelihood.
- MR conditional strength: the old conditional-F accessor is deprecated and
  throws because it only held marginal F values. True conditional strength needs
  cross-exposure covariance; use the explicitly named marginal accessor meanwhile.
- Selective inference: unknown-noise or response-selected penalties, dependent
  sample splits, and non-Gaussian selection require different inferential methods.
- Quantile/nonparametric regression: quantile inferential covariance,
  mixed-type/multidimensional kernels, and automatic inference-valid bandwidth
  selection remain distinct extensions; exact nonsmoothed quantile fitting is implemented.
- Kernel set tests: ill-conditioned positive-mixture spectra can exceed the
  bounded gamma-series calculation and now fail explicitly. Analytic SKAT-O
  remains moment-matched; no exact finite-sample calibration is implied.
- Numerical identification: extreme original-coordinate deficient SVD fits
  and inconclusive sparse covariance-rank checks reject rather than returning
  a numerically unsupported estimator or variance decomposition.
