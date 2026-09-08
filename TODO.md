# Development inventory

Last reviewed: 2026-09-07.

This inventory keeps unfinished decisions at the top. Completed work is archived
at the bottom so it remains searchable without obscuring the active backlog.

Priority meanings:

- **Major**: removes a material scale limit or adds a cross-cutting analysis
  path.
- **Medium**: adds an important inference, workflow, or documentation feature
  without changing the core architecture.
- **Candidate**: a documented boundary that still needs a product decision.

## Active committed work

No major or medium items are currently in progress.

## Requested future analyses

No requested analysis items remain in this inventory. The six regression
families requested in the prior review are implemented below with explicit,
deterministic contracts and corresponding tests and vignettes.

## Candidates requiring prioritization

These are real documented boundaries that remain outside the bounded extension
contracts below:

- full latent-variable RAM optimization with ordinal DWLS/polychoric
  likelihoods and multi-group SEM;
- publication-bias models beyond the deterministic PET/PEESE and labelled
  median-reflection diagnostic;
- high-dimensional unstructured sparse covariance optimization beyond the
  one-block coordinate-refit path;
- state-space diffuse likelihoods with arbitrary seasonal missingness beyond
  the exposed differenced and stationary observed-pattern paths;
- generalized MR estimators requiring additional summary-data fields beyond
  LD-correlated IVW/Egger and overlap-aware IVW;
- multi-dimensional adaptive quadrature for GLMMs; and
- selective-inference polyhedral truncation after the deterministic split.

Sources: [SEM](docs/vignettes/sem.md),
[meta-analysis](docs/vignettes/meta-analysis.md),
[mixed-model roadmap](docs/lme4-pedigreemm-roadmap.md),
[time series](docs/vignettes/time-series.md),
[end-to-end MR](docs/vignettes/mr-end-to-end.md),
[GLMM/pedigree](docs/vignettes/pedigree-and-glmm.md), and
[linear/penalized models](docs/vignettes/linear-models-and-glms.md).

## Boundaries that are not automatically TODO items

Some vignette warnings state a statistical contract rather than a missing
feature. Do not turn these into implementation work without a separate design
decision. Examples include PQL not being an exact marginal GLMM likelihood,
ordinary post-selection OLS p-values not being selection-adjusted, MR
sensitivity analyses not proving the exclusion restriction, and the xWAS
benchmark not measuring an end-to-end TwoSampleMR workflow.

## Completed work

### Advanced analysis extensions — completed

Completed 2026-09-07. The seven previously documented candidate areas now have
bounded, tested implementations. SEM adds principal-factor latent measurement,
mean/intercept summaries, marginal ordinal thresholds, Gaussian FIML missing
patterns, score-based sandwich/cluster corrections, modification-index
diagnostics, and indirect-effect delta-method inference. Meta-analysis adds
full-covariance GLS moderator regression plus PET, PEESE, and labelled
median-reflection trim-and-fill diagnostics. Mixed models add sparse
unstructured covariance-shape coordinate refits and a formula-native entry
point. Time series add explicit diffuse-state ARIMA bookkeeping and exact
stationary observed-pattern handling for missing values. MR adds generalized
and overlap-aware CLI methods with LD-based conditional inference. GLMM adds a
fixed 10-node Gaussian-quadrature random-intercept path. Penalized regression
adds deterministic held-out data-splitting inference after LASSO/elastic-net
selection.

The remaining boundaries are intentionally listed above: these implementations
are not claims of full lavaan, metafor, lme4, forecast, TwoSampleMR, AGQ, or
polyhedral-selective-inference parity.

Source: [advanced extensions vignette](docs/vignettes/advanced-extensions.md).

### Regression families — completed

Completed 2026-09-07. The regression package now provides shared-design
Gaussian multivariate OLS, baseline-category multinomial logistic regression,
smoothed pinball-loss quantile regression, Gaussian-kernel nonparametric
regression, deterministic span-selecting supersmoothing, and a partially
linear semiparametric model. These paths are in-memory Java APIs and do not
use bootstrap or Monte Carlo sampling. `RegressionFamiliesBenchmark` reports
repeatable JVM throughput for all six families.

The R relationship is intentionally documented rather than overstated:
multivariate OLS is the direct Gaussian reference case, while multinomial,
supersmoother, quantile, nonparametric, and semiparametric implementations
are compatible contracts rather than byte-for-byte reimplementations of
`nnet`, `stats::supsmu`, `quantreg`, `np`/`npreg`, or a single semiparametric
R package. Direct R timing was not run on this host because `Rscript` is not
installed; the benchmark output is therefore a Java baseline, not a claimed
R speedup.

Source: [regression families vignette](docs/vignettes/regression-families.md).

### Nonlinear fixed and mixed effects — completed

Completed 2026-09-07. `NonlinearFixedModel` fits Gaussian nonlinear means with
analytic gradients using damped Gauss–Newton iterations. `NonlinearMixedModel`
linearizes the same mean function while reusing prepared sparse REML equations
for ordinary random effects, and `fitPedigree` adds sparse numerator-relationship
precision terms with optional ordinary effects. Results expose nonlinear
parameters, fitted values, residuals, objective, convergence metadata, and the
final sparse linearized mixed-model result.

The implemented contract is deliberately explicit: this is an additive
Gaussian nonlinear mixed-effects path, not a non-Gaussian Laplace/AGQ solver or
an exact `nlme` likelihood implementation. Analytic gradients are required,
and callers should validate parameterizations against independent reference
fixtures.

Source: [nonlinear fixed and mixed-effects vignette](docs/vignettes/nonlinear-models.md).

### Meta-analysis, sparse correlated mixed models, and MR workflow extensions — completed

Completed 2026-09-07. Meta-analysis now has first-order effect-size
construction for mean differences, Hedges g, log odds/risk ratios, and Fisher-z
correlations; generalized inverse-covariance pooling for correlated study
effects; cluster-robust meta-regression sandwich inference; and deterministic
Egger/rank publication-bias diagnostics. Mixed models now expose sparse
grouped correlated blocks with caller-supplied within-group covariance shapes,
a sparse Satterthwaite contraction utility, and a generic profile-likelihood
interval solver. The existing `lme4`-style parser continues to provide
intercepts, slopes, `||`, nested grouping, and correlated single-bar blocks.
MR now includes `mr-estimate` for the core independent estimators, dependency-
free native SVG scatter plots, and an explicit secondary-signal clumper for
conditional LD selection.

These additions are intentionally bounded: correlated meta-analysis currently
fits an intercept-only GLS model; sparse correlated mixed blocks estimate one
scale per supplied covariance shape; the finite-DF utility consumes derivatives
from a caller's sparse factorization; and conditional MR selection requires an
LD correlation matrix rather than claiming genotype-level conditional p-values.

Completion evidence: focused JUnit coverage for effect-size construction,
identity-covariance parity, robust/publication-bias diagnostics, sparse
correlated fitting, profile crossings, native SVG output, and secondary-signal
LD selection.

### Frequentist linear mediation analysis — completed

Completed 2026-09-07. `MediationAnalysis` fits the mediator, outcome, and
total-effect Gaussian OLS models and reports the `a`, `b`, indirect, direct,
and total effects. Sparse REML mediation also supports ordinary random-effect
terms and additive pedigree precision terms, reusing prepared sparse
structure across all three fits. Indirect-effect inference uses the analytic
Sobel/product-of-coefficients delta method with a normal confidence interval;
no bootstrap, Monte Carlo, or posterior sampling is used. Common complete-row
selection keeps all three component fits on the same analysis sample.

Source: [mediation vignette](vignettes/mediation.md).

### Frequentist zero-inflated mixed and pedigree models — completed

Completed 2026-09-05. `SparseZeroInflatedMixedModel` now fits ZIP and NB2-ZINB
models with sparse random effects in both the count and structural-zero
predictors. For zero observations the random-mode system includes the full
observed count, zero, and cross curvature. Independent terms accept identity or
caller-supplied precision matrices; `CorrelatedZeroInflatedRandomEffect`
estimates two variances and a guarded correlation with joint precision
`inverse(G) kron inverse(A)`.

BOBYQA remains derivative-free and sees only the low-dimensional outer
parameters. Sparse damped Newton solves the random mode. The optional
`fitWithInference` path computes a post-fit numerical observed Hessian of the
marginal Laplace objective, including nuisance covariance parameters, and
reports fixed/count/zero/dispersion covariance and standard errors. Prepared
fits reuse one symbolic sparse analysis and one numerical factor per worker;
profile-likelihood grids and deterministic conditional parametric bootstrap
reuse that prepared structure. Results expose boundary/singularity warnings.

Completion evidence:

- analytic two-predictor scores and observed Hessian blocks match finite
  differences for ZIP and ZINB zero and positive observations;
- ordinary and pedigree effects work independently in either process, and a
  simulation-recovery test resolves structural-zero group ordering;
- compiled count/zero/size formulas avoid repeated parsing and design builds;
- grouped estimates, marginal log likelihoods, variance components, and ZIP
  standard errors pass checked-in `glmmTMB` 1.1.14 fixtures;
- automatic bounded BFGS uses parallel full-Laplace numerical gradients for
  small modes while large sparse modes and failed BFGS fits retain BOBYQA;
- a separately compiled TMB 1.9.25 sparse-GMRF template gates correlated
  pedigree fixed effects, variances, correlation, and marginal likelihood;
- the optimized reproducible 500-member, 1,500-observation benchmark used
  1,000 random coefficients, 4,900 sparse equation/factor nonzeros, converged
  in 2.161532 seconds (previously 2.576353), and sampled a 312,168,592-byte
  peak-heap delta on the documented development host.

Sources: [pedigree and GLMM vignette](docs/vignettes/pedigree-and-glmm.md),
[`SparseZeroInflatedMixedModel`](src/main/java/org/jlinalg/distributional/SparseZeroInflatedMixedModel.java),
and the [R/TMB generators](src/test/resources/r-reference).

### Sparse Cox mixed and pedigree solves — completed

Completed 2026-09-05. Generic sparse precision, repeated observations,
prepared symbolic factorization, sparse diagnostics, dense/sparse parity, and a
deterministic `coxme` performance gate are implemented. The 800-row reference
scan was 10.03x faster with eight JLinAlg workers than single-threaded R.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Scalable pedigree prediction uncertainty — completed

Completed 2026-09-05. Sparse LMM diagonal PEV/reliability extraction uses
batched inverse-column solves, supports selected members, agrees with the dense
reference, and crosses the former 256-coefficient cutoff. Full dense covariance
remains intentionally quadratic rather than part of this contract.

Source: [pedigree vignette](docs/vignettes/pedigree-and-glmm.md).

### Prepared Cox score-scan pipeline — completed

Completed 2026-09-05. One null model and risk-set plan support deterministic
bounded, file-backed scans with model-based, cluster-robust, or supplied
correlation variance. The checked-in fixture measured 28,748 predictors/second.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### xWAS all-pairs output and BH/FDR correction — completed

The xWAS pipeline can stream all analyzable pairs, applies disk-backed BH over
the complete prespecified finite-p family, and records exclusions, failures,
test counts, and adjusted values explicitly.

Source: [xWAS pipeline](docs/vignettes/xwas-mr-pipeline.md).

### xWAS hit hand-off and robust follow-up — completed

Completed 2026-09-05. Retained hits preserve harmonized instruments and feed
RAPS, contamination-mixture, and PRESSO-style follow-up with deterministic
ordering, warnings, partial results, and CLI persistence.

Source: [xWAS follow-up](docs/vignettes/xwas-mr-pipeline.md#follow-up-analyses-for-hits).

### Dedicated colocalization workflow and MR documentation — completed

Completed 2026-09-05. Markdown and website vignettes cover alignment, priors,
H0-H4, shared-variant posterior, diagnostics, and MR/xWAS hand-off, with
navigation and stale statements repaired.

Source: [colocalization vignette](docs/vignettes/colocalization.md).

### Cox inference and diagnostics — completed

Completed 2026-09-05. Cluster sandwich covariance, residual/influence exports,
term/global proportional-hazards tests, and fixed/profiled gamma frailty are
available and independently gated.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Start-stop Cox performance path — completed

Completed 2026-09-05. An indexed ascending sweep maintains delayed-entry risk
moments while preserving Breslow/Efron behavior and the reference path.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Incremental omics result sink — completed

Confirmed 2026-09-05. Association and Cox scans stream deterministic estimates,
failures, and accounting rather than retaining large result grids in memory.

Source: [GWAS/TWAS pipeline](docs/gwas-twas-pipeline.md).

### Zero-inflated outer-point cache optimization — completed

Completed 2026-09-05. Worker-local caches remove repeated fixed-predictor,
NB2 gamma/size, and eligible zero-probability calculations from inner mode
solves. Full Salamanders fits improved another 1.41x (ZIP) and 1.83x (ZINB),
with identical reported estimates and standard errors. Added reproducible
Java/R runners, raw paired timings, independent likelihood reconstruction,
and parallel/cache-refresh tests. The large pedigree case showed no speedup.

Source: [benchmark and validation report](docs/zero-inflated-cache-performance.md).
