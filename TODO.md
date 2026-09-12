# Development TODO

Last reviewed: 2026-09-12.

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

## Rare-variant meta-analysis — further scope

Implemented engineering, calibration, conditional/VT/heterogeneous tests, and
leave-out diagnostics are documented in the
[rare-variant tutorial](docs/rare-variant-meta-analysis.md) and
[validation report](docs/rare-meta-validation.md).

The Raremetal2 format/method assessment and quantitative model-metadata contract
are completed in the [trait-model assessment](docs/raremetal2-trait-models.md).
The remaining implementations are separate:

- [ ] **Allele-aware multiallelic summaries.** Implement producer-specific
  compressed covariance decoding and allele-aware union/group/conditioning
  keys, preserving both covariance axes and global row ordering. Validate
  same-site cross-ALT covariance, missing coverage and indexed/parallel equality.
  Old position-only covariance cannot resolve multiallelic groups.
- [ ] **Unbalanced-study exact method.** Specify the pooled estimand and required
  phenotype moments, residualization and covariance scaling before implementing
  Raremetal2 `--useExact`. Validate against independent participant-level fits
  and cross-products; this is separate from binary rare-case calibration.
- [ ] **Additional rare-meta trait models and calibration.** Implement a richer
  versioned model protocol and validated binary rare-case/joint group tails,
  followed by model-specific Poisson/survival support. Reuse the existing
  conditional-GWAS score infrastructure where appropriate; normal-only exports
  and cohort p-values do not supply calibrated rare-event group inference.

## Validation for new methods

Use checked-in fixtures from the underlying R estimators, independent
likelihood/derivative/covariance checks, and reproducible workload-specific
benchmarks. Zelig often wraps other packages and must not be the sole numerical
reference. Document estimator assumptions, convergence failures, and supported
scope; scenario contrasts require additional assumptions for causal interpretation.

## Conditional GWAS — further scope

The opt-in Gaussian/logistic/Poisson/Cox aggregate export and cohort-side
conditioning refits are implemented; see the
[summary schema and validation](docs/conditional-gwas-summary.md).

- [ ] **Model-specific tail calibration.** Add validated binary rare-case and
  survival calibration, with the required score cumulants/risk-set information.
  Current exports explicitly provide normal tails and no calibrated p-value.
- [ ] **Summary-only score import and conditioning.** Add an importer for the
  new schema with allele alignment, null-model compatibility, full matrix
  coverage and rank checks. Missing cross-block covariance must remain unknown.
  Keep local Schur-complement inference separate from nonlinear null refits;
  the existing `mr-estimate` Gaussian summary interface is unchanged.
- [ ] **Additional cohort score models.** Add dedicated mixed, cluster-robust,
  stratified/repeated-subject Cox, grouped-binomial and other outcome exporters.
  Define cross-cohort overlap handling before pooling dependent study scores.

## Remaining estimator extensions

- [ ] **SEM.** DWLS/WLSMV, mixed continuous/ordinal responses, general MAR
  ordinal missingness, multigroup invariance, ordinal modification indices,
  and FIML robust scaled fit statistics.
- [ ] **Time series.** General regression terms and historical state smoothing.
- [ ] **Nonparametric inference.** Mixed-type/multidimensional kernels and
  automatic inference-valid bandwidth selection; dependent observations and
  response mass points require separate quantile inference methods.
- [ ] **Non-tensor random-effect integration.** Add a non-tensor integration
  method before expanding quadrature to higher-dimensional structures.

## Interpretation and numerical limits

These boundaries are explicit; this release does not claim complete parity
with every feature of lavaan, metafor, lme4, or the other R packages.

- The legacy MR conditional-GWAS interface assumes its documented Gaussian
  score model. Cohort model-specific scores/covariance support local inference;
  they do not reconstruct a complete nonlinear likelihood away from the null.
  New nonlinear conditioning sets require cohort-side refits or a richer
  likelihood evaluation protocol.
- Multidimensional quadrature is deliberately tensor-node bounded and supports
  binomial-logit and Poisson-log likelihoods. Higher-dimensional structures
  should use sparse Laplace/PQL unless a non-tensor integration method is added.
- Zero-inflated deterministic multi-start certification is relative to the
  configured starts, not a mathematical global optimum over the parameter
  continuum. The random-effect integral remains a single-mode Laplace
  approximation rather than a sum over separated conditional modes.
- Nonlinear models: floating-point quantization under very large response
  offsets can prevent a score certificate even when fitted SSE is near optimal.
  Such fits return nonconvergence. Mixed effects are additive Gaussian effects,
  not a general nonlinear random-parameter likelihood.
- Selective inference: unknown-noise or response-selected penalties, dependent
  sample splits, and non-Gaussian selection require different inferential methods.
- Exact quantile covariance uses supplied conditional densities or an iid
  residual kernel with caller-set bandwidth. This does not establish automatic
  bandwidth validity or inference for dependent rows/response mass points.
- Kernel set tests: the saddlepoint positive-mixture fallback is approximate,
  as is opt-in analytic SKAT-O. Default parametric SKAT-O simulation is exact
  for its Gaussian score-null model up to Monte Carlo error, not an exact
  finite-sample phenotype calibration.
- Numerical identification: extreme original-coordinate deficient SVD fits
  and inconclusive sparse covariance-rank checks reject rather than returning
  a numerically unsupported estimator or variance decomposition.
