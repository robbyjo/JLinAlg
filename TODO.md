# Development TODO

Last reviewed: 2026-09-19.

This inventory tracks open work. Completed implementation and audit details are
recorded in the [advanced-method validation report](docs/advanced-validation.md),
the [v0.3.0 remaining-code audit](docs/release-0.3.0-audit.md), the
[release notes](RELEASE_NOTES.md), and the
[regression inference validation report](docs/regression-inference-validation.md).

## xWAS analysis additions

Implemented differential, regional, multiplicity and imputation workflows are
documented in [modern inference validation](docs/inference-workflows-validation.md).

The first four follow-up workflows (unpartitioned LDSC, genetically predicted
TWAS/PWAS, a single genetic factor with conditional SNP tests, and prediction
score training/application/evaluation) are implemented in the source build.
See [xWAS workflow validation and scope](docs/xwas-followup-validation.md).

- [ ] **Frozen confounder projection and count-aware adjustment.** Add
  fold-owned fit/freeze/apply artifacts for prediction, plus independently
  validated RUV and ComBat-Seq contracts. Do not refit preprocessing on held-
  out folds or apply Gaussian ComBat directly to raw counts.

Further extensions of the new workflows remain explicit scope boundaries:
partitioned/two-step/liability-scale LDSC and native summary munging; native
PredictDB/FUSION model adapters and rank-truncated multi-tissue inference;
multiple genomic factors, robust DWLS and propagation of measurement-model
uncertainty into SNP effects; LD-aware Bayesian PRS, logistic score training,
grouped CV and absolute-risk calibration. Current methods do not claim full
parity with LDSC, MetaXcan, GenomicSEM, or PRS software suites.

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

The opt-in Gaussian/logistic/probit/Poisson/Cox aggregate export, cohort-side
conditioning refits, and summary-only importer are implemented; see the
[summary schema and validation](docs/conditional-gwas-summary.md).

- [ ] **Model-specific tail calibration.** Add validated binary rare-case and
  survival calibration, with the required score cumulants/risk-set information.
  Current exports explicitly provide normal tails and no calibrated p-value.
- [ ] **Additional cohort score models.** Add dedicated mixed, cluster-robust,
  stratified/repeated-subject Cox, grouped-binomial and other outcome exporters.
  Define cross-cohort overlap handling before pooling dependent study scores.

## Remaining estimator extensions

Implemented methods and their tested scope are documented in the
[estimator extension report](docs/estimator-extensions.md). Exact diffuse
historical smoothing now scales to long series, and ARIMA regression exposes
both conditional forecasts and forecasts with joint regression/dynamic
parameter uncertainty.

- [ ] **SEM scalability and additional test variants.** Extend joint ordinal
  integration beyond four observed ordinal responses per row; add multigroup
  WLSMV, DWLS-specific modification indices, robust nested-model comparisons,
  and robust CFI/TLI. Current multigroup and ordinal modification tests use the
  joint likelihood, with explicit equality constraints and MAR marginalization.
- [ ] **Discrete conditional quantile inference.** Extend beyond the iid
  marginal order-statistic interval to regression coefficients for responses
  with mass points. Cluster/HAC density-based regression inference still assumes
  continuous responses. Automatic bandwidths do not promise finite-sample or
  simultaneous coverage.

## Interpretation and numerical limits

These boundaries are explicit; this release does not claim complete parity
with every feature of lavaan, metafor, lme4, or the other R packages.

- The legacy MR conditional-GWAS interface assumes its documented Gaussian
  score model. Cohort model-specific scores/covariance support local inference;
  they do not reconstruct a complete nonlinear likelihood away from the null.
  New nonlinear conditioning sets require cohort-side refits or a richer
  likelihood evaluation protocol.
- Multidimensional integration supports binomial-logit and Poisson-log models.
  Tensor Hermite remains node bounded. The optional defensive randomized-QMC
  method uses empirical error estimates and bounded dense precision operations;
  it is not certified integration or a replacement for large sparse methods.
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
- Exact quantile regression has iid residual-density, supplied-density,
  cluster, and HAC inference paths with distinct regularity assumptions.
  Response mass points use a separate marginal interval API; conditional
  quantile-regression inference at mass points remains open.
- Kernel set tests: the saddlepoint positive-mixture fallback is approximate,
  as is opt-in analytic SKAT-O. Default parametric SKAT-O simulation is exact
  for its Gaussian score-null model up to Monte Carlo error, not an exact
  finite-sample phenotype calibration.
- Numerical identification: extreme original-coordinate deficient SVD fits
  and inconclusive sparse covariance-rank checks reject rather than returning
  a numerically unsupported estimator or variance decomposition.
