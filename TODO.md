# Development inventory

Last reviewed: 2026-09-05.

This inventory collects explicit future-work statements from the documentation
and website vignettes, plus gaps found while reviewing the implemented APIs. It
tracks product priority separately from estimated effort: a small change can be
important, and a large change can remain deferred.

Priority meanings:

- **Major**: removes a material scale limit or adds a cross-cutting analysis
  path.
- **Medium**: adds an important inference, workflow, or documentation feature
  without changing the core architecture.
- **Candidate**: documented boundary that needs a product decision before it is
  promoted to committed work.

## Major priority

### Sparse Cox mixed and pedigree solves — completed

Completed 2026-09-05. `SparseCoxMixedModel` now accepts a generic
unit-incidence term plus caller-supplied `SparsePrecisionMatrix`, supports
repeated observations per coefficient, and exposes prepared symbolic
factorization. `CoxPedigreeFrailty.fitSparse` accepts either an existing
`Pedigree` or an end-to-end sparse `PedigreeRandomEffectTerm`.
Results identify the dense or sparse solver, report sparse equation/factor
nonzeros and lower-bound singular fits, and retain convergence/backend
provenance. The TOPMed Cox benchmark now records sampled peak heap and sparse
matrix diagnostics. Dense/sparse tests verify conditional-mode parity,
including repeated pedigree IDs. A deterministic 800-row repeated-measures
fixture also gates coefficients, standard errors, convergence, and elapsed time
against R `coxme`; the 2026-09-05 run passed and was 10.03x faster with eight
JLinAlg scan workers than the single-threaded R scan.

Completion criteria met:

- the sparse path does not materialize a dense random-effect covariance or
  information matrix;
- dense and sparse conditional estimates agree on numerical fixtures;
- singularity, convergence, solver choice, sparse structure, and backend
  behavior are reported explicitly; and
- the checked-in benchmark records time, sampled peak heap, numerical parity,
  and an executable R comparison gate.

The completed kernel deliberately remains an approximation: it retains the
exact penalized score but uses the diagonal of the profiled random-effect
information with the full sparse precision. Exact sparse Laplace
determinants/covariance, multiple sparse terms, ties/strata/start-stop support,
a sparse-GRM convenience
facade, and private-cohort large-pedigree reruns are follow-on enhancements,
not blockers for this completed item.

Sources: [Cox vignette](docs/vignettes/cox-survival.md),
[website survival vignette](site/vignettes/survival.html).

### Scalable pedigree prediction uncertainty — completed

Completed 2026-09-05. Sparse LMM diagonal PEV extraction now solves inverse
columns in batches of 32 rather than disabling uncertainty above 256 random
coefficients. `SparsePedigreeRemlResult` exposes all PEVs/reliabilities and
selected-member lookup. Small-pedigree values match the dense reference; a
300-member regression test crosses the former cutoff. The checked-in benchmark
measured 500 finite PEVs from 1,000 observations in 0.076384 seconds on the
documented development host. A full covariance is deliberately not promised:
its output alone is quadratic.

Provide prediction-error variance and reliability from the sparse combined
pedigree model. The API and documentation should use population-neutral terms
such as individuals or pedigree members; animal breeding is one use case, not
the only one.

A complete conditional covariance is inherently quadratic in the number of
predicted individuals. The scalable target should therefore be diagonal PEV and
reliability for all requested individuals, with optional selected covariance
blocks, rather than promising that an arbitrarily large full dense matrix is
cheap.

Done when:

- diagonal PEV/reliability can be extracted without forming the full inverse;
- callers can request selected individuals or covariance blocks;
- results agree with the dense reference path on small pedigrees; and
- memory scaling is documented and benchmarked.

Sources: [pedigree vignette](docs/vignettes/pedigree-and-glmm.md),
[website pedigree vignette](site/vignettes/pedigree.html), and the
[compatibility roadmap](docs/lme4-pedigreemm-roadmap.md).

### Prepared Cox score-scan pipeline — completed

Completed 2026-09-05. `FastCoxAssociation` fits one null model and reuses the
risk-set plan for deterministic bounded blocks with model-based,
cluster-robust, or caller-correlation variance. File-backed omics scans stream
through `scanPredictorsCoxTo` and the existing sink/failure accounting. The
checked-in 2,000-by-512 benchmark measured 28,748 predictors/second versus
2,202 predictors/second for Java full Cox refits and 1,829 predictors/second
for conventional R `survival::coxph.fit` refits on the same host. The prepared
path was 15.72x faster than R for this representative throughput fixture.

Add a prepared null-model score scan for many genetic or molecular predictors,
including file-backed blocks and the existing cohort/QC conventions. This is
the survival-analysis analogue of the association pipelines and should reuse a
fixed Cox risk-set plan rather than refitting the null model for every marker.

Done when:

- right-censored score tests reuse one fitted null model and risk-set plan;
- block streaming, deterministic output, failures, and provenance match the
  association pipeline conventions;
- robust and relatedness-aware variance can be selected when those engines are
  available; and
- numerical fixtures and a representative high-throughput benchmark are
  checked in.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

## Medium priority

### Complete frequentist zero-inflated mixed and pedigree models — in progress

The first vertical slice landed 2026-09-05. `SparseZeroInflatedMixedModel`
fits zero-inflated Poisson and NB2 models with separate fixed count,
structural-zero, and NB2-size designs. One or more ordinary or pedigree random
effects may enter the count process. BOBYQA optimizes the non-random parameters;
the high-dimensional random-effect mode uses analytic observed curvature,
sparse damped Newton steps, and the same Hessian for the first-order Laplace
determinant. Pedigree fits consume Henderson `A^-1` directly and retain
unphenotyped ancestors.

The following work remains:

1. admit random effects in the structural-zero predictor and assemble the full
   count/zero cross-Hessian for zero observations;
2. support independent and correlated two-process pedigree effects, with the
   latter using the Kronecker precision `inverse(G) kron inverse(A)` and guarded
   correlation transforms;
3. add the observed numerical Hessian of the marginal Laplace objective,
   fixed-effect/dispersion covariance, profile likelihood, and parametric
   bootstrap for boundary-prone variance and zero-inflation parameters;
4. add a prepared repeated-fit API with one symbolic sparse factorization and
   one numeric factor per worker;
5. gate grouped ZIP/ZINB values against `glmmTMB` and pedigree values against an
   independent frequentist TMB sparse-GMRF template; and
6. add formula compilation, simulation/recovery tests, singularity diagnostics,
   and representative sparse pedigree performance/heap benchmarks.

Done when both distributional processes can carry sparse random effects,
correlated pedigree structure has an explicit estimable covariance contract,
reported marginal inference is independently validated, and prepared large-
pedigree fitting is demonstrably sparse and reproducible.

Sources: [pedigree and GLMM vignette](docs/vignettes/pedigree-and-glmm.md) and
[`SparseZeroInflatedMixedModel`](src/main/java/org/jlinalg/distributional/SparseZeroInflatedMixedModel.java).

### xWAS all-pairs output and BH/FDR correction — completed

The xWAS MR pipeline retains threshold-passing pairs for full diagnostics and
can now stream every successfully screened pair to an optional all-pairs table.
The CLI appends a BH-adjusted value across the complete analyzable family.

The implementation reuses the CLI's disk-backed `ExternalBh` path, so the
bounded scan does not retain the complete exposure-by-outcome grid in memory.

Completed behavior:

- BH uses every finite screening p-value in the prespecified exposure-by-outcome
  family, without a nominal-p prefilter;
- excluded, insufficient-instrument, and failed pairs have explicit semantics
  and are not silently counted as tests;
- the CLI records the number of adjusted tests and emits `fdr_bh`; and
- adjusted values match an independent R/JDistlib reference including ties,
  zero, one, and non-finite results.

Sources: [xWAS pipeline vignette](docs/vignettes/xwas-mr-pipeline.md) and
[`ExternalBh`](src/main/java/org/jlinalg/cli/ExternalBh.java).

### xWAS hit hand-off and robust follow-up — completed

Completed 2026-09-05. Only retained hits keep their immutable harmonized
instrument list. `XwasMrBatchResult.followUp` runs RAPS, contamination mixture,
and PRESSO-style analysis in deterministic hit order while preserving method
warnings and partial results. `mr-xwas --follow-up-output` persists estimates,
convergence, likelihood/probability, global p-values, distortion, outliers,
and warnings without rerunning the exposure-outcome grid.

Make retained hits easy to pass into RAPS, contamination-mixture, and
PRESSO-style analyses without redoing the full scan. Prefer a bounded design:
retain or reconstruct harmonized instruments only for hits rather than storing
them for every pair. Decide whether follow-up runs automatically or through a
second explicit CLI stage, and persist convergence, outliers, warnings, and
method-specific p-values.

Source: [xWAS hit follow-up](docs/vignettes/xwas-mr-pipeline.md#follow-up-analyses-for-hits).

### Dedicated colocalization vignette and MR documentation repair — completed

Completed 2026-09-05. The standalone Markdown and synchronized website
vignettes cover alignment, priors, overlap trimming, H0-H4, conditional shared
variant posterior, diagnostics, and MR/xWAS hand-off. Both vignette indexes
link the new page and the end-to-end MR guide now points to it.

`ColocSusie` is implemented, tested against `coloc::coloc.susie`, benchmarked,
and briefly demonstrated inside the combined SuSiE/SEM vignette. What is
missing is a dedicated colocalization workflow. Add a vignette covering input
alignment, priors, posterior overlap/trimming, H0-H4 interpretation,
variant-level shared posterior, diagnostics, and an MR/xWAS follow-up example.

At the same time, remove the stale end-to-end MR statement that JLinAlg does not
implement colocalization, add the new page to both vignette indexes/navigation,
and keep the generated website synchronized with the Markdown source.

Sources: [current compact example](docs/vignettes/susie-and-sem.md),
[stale MR statement](docs/vignettes/mr-end-to-end.md), and
[colocalization implementation](src/main/java/org/jlinalg/coloc/ColocSusie.java).

### Cox inference and diagnostics — completed

Completed 2026-09-05 as separate APIs: `CoxDiagnostics` provides cluster-
robust covariance plus martingale, deviance, score, dfbeta, and Schoenfeld
exports; `CoxProportionalHazardsTest` provides term/global log-time score tests;
and `CoxGammaFrailty` provides multiplicative shared gamma frailty with fixed
or Laplace-profiled variance. Fixed Cox coefficients, cluster sandwich SEs,
and fixed-theta gamma estimates are gated against R `survival` fixtures.

Add the following as separable deliverables so they can land independently:

1. cluster-robust sandwich covariance, including recurrent-event use;
2. proportional-hazards diagnostics based on Schoenfeld/score residuals and a
   documented global/test-per-term interface;
3. influence and residual exports suitable for diagnostic plots; and
4. gamma frailty, with the approximation and integration method explicit.

These should be validated independently rather than hidden behind one broad
"Cox diagnostics" completion flag.

Sources: [Cox vignette](docs/vignettes/cox-survival.md),
[website survival boundary](site/vignettes/survival.html).

### Start-stop Cox performance path — completed

Completed 2026-09-05. `CoxCountingProcessPlan` indexes starts, stops, and event
groups once per stratum. An ascending sweep maintains risk moments while rows
enter and leave, preserving Efron/Breslow behavior and the independent
statsmodels delayed-entry fixture. The checked-in benchmark measured a
2,000-row start-stop fit in 0.004636 seconds on the development host.

Start-stop observations currently use the general risk-set reference path.
Profile representative recurrent-event and time-dependent-covariate workloads,
then add an indexed/sweep implementation if the benchmark confirms a material
bottleneck. Preserve the reference path for parity tests.

Source: [Cox vignette](docs/vignettes/cox-survival.md).

### Incremental omics result sink — completed

Completed before this review and confirmed 2026-09-05.
`OmicsAssociationSink`, `OmicsAssociationSummary`, `scanPredictorsTo`, and
`scanPredictorsGlmTo` already provide deterministic incremental estimates,
failures, and accounting; the stale guide was repaired. The new Cox streaming
entry point uses the same contract.

Complete the planned incremental sink for omics response scans so large
feature-by-predictor results need not be retained in memory. Match the existing
association sink's accounting, deterministic ordering, and failure contract.

Source: [GWAS/TWAS pipeline guide](docs/gwas-twas-pipeline.md).

## Candidates requiring prioritization

These are real documented boundaries, but the vignettes do not establish that
they are committed major- or medium-priority work:

- **SEM:** latent measurement variables, mean/intercept structures, ordinal
  thresholds, robust/clustered corrections, modification indices,
  indirect-effect delta-method inference, and FIML missing-data patterns.
- **Meta-analysis:** cluster-robust variance, multilevel/correlated effects,
  publication-bias diagnostics, and effect-size construction.
- **Mixed models:** sparse correlated-block likelihoods, scalable sparse
  finite-DF calculations, broader `lme4` formula parity, and profile-likelihood
  intervals.
- **Time series:** diffuse exact likelihood for integrated models and a sparse
  missing-series path.
- **MR workflow:** a general estimator CLI, native plot rendering, and
  conditional/secondary-signal clumping.
- **GLMM:** exact or adaptive-quadrature alternatives where PQL is inadequate.
- **Penalized regression:** selection-aware inference after LASSO/elastic-net
  selection.

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
