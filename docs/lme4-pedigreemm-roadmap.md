# lme4 and pedigreemm functional roadmap

JLinAlg preserves the current matrix-first OLS, GLM, dense REML, pedigree REML,
and PQL APIs. The compatibility layer is additive: callers that already supply
covariance bases continue to receive the same numerical model.

## Implemented foundation

- Common coefficient association results: beta, SE, Wald t/z statistic,
  denominator DF, and two-sided p-value.
- Fast default mixed-model DF `N - rank(X) - 1`.
- Optional coefficient-specific, delta-method Satterthwaite DF based on the
  REML variance-parameter information and analytic `vcov(beta)` derivatives.
- Optional Kenward-Roger covariance adjustment and coefficient-specific DF for
  linear covariance-component REML models.
- Random-intercept, independent random-slope, and general `Z Z'` covariance
  constructors.
- A general multi-term Gaussian REML facade with crossed/nested independent
  terms, conditional random-effect modes, PEVs, fitted values, and residuals.
- Sparse-CSR grouped random-intercept and random-slope design storage, with a
  dense compatibility view for existing callers.
- Pedigree Gaussian REML and a pedigree binomial/Poisson PQL facade.
- Direct sparse-CSR construction of pedigree `A^-1` without dense numerical
  inversion.
- JDistlib backend selection and provenance for all numerical fitting paths.
- Profile ML as an alternative to REML, joint linear-contrast tests,
  singular-boundary diagnostics, and diagonal residual-weight bases.
- Compiled fixed and mixed formulas with treatment/sum contrasts, interactions,
  offsets/weights, and sparse independent random intercept/slope terms.
- Henderson prediction equations consuming sparse grouped `Z`; pedigree
  BLUP/PEV/reliability now consumes sparse `A^-1` directly.
- Batched P3D/EMMAX-style null-model reuse for GWAS/TWAS marker scans.
- Related-sample Burden, SKAT, and SKAT-O score tests reuse the same retained
  REML projection rather than refitting a mixed model for every set.
- Sparse-equation ML/REML for independent grouped terms with reusable
  minimum-degree sparse Cholesky and no observation-scale covariance matrix.
- Analytic-gradient Gaussian nonlinear fixed and mixed effects, with prepared
  sparse ordinary and pedigree random-effect linearizations.
- Sparse pedigree variance estimation using `A^-1` directly.
- Cholesky-parameterized unstructured correlated random blocks on the dense
  reference likelihood path.
- Formula `||` independent terms, nested grouping shorthand, correlated
  single-bar blocks, and `fixef`/`ranef`/`VarCorr`-style accessors.
- ML-only nested likelihood-ratio model comparison.
- Retained-design response refits with variance-component warm starts.
- Marginal and conditional Gaussian response simulation for ordinary LMMs and
  pedigree animal models, with deterministic seeds.
- Accelerator-friendly parametric bootstrap for fixed effects and variance
  components, including percentile intervals, empirical bias/SE, bounded CPU
  parallelism, warm refits, and structured convergence failures.
- Marginal and conditional prediction on new data, including an explicit
  zero-mode policy for unseen grouping levels.
- Multiple named pedigree precision terms and ordinary independent random
  terms in the same sparse REML model, including unphenotyped ancestors.

## High-performance mixed-model core

The sparse engine represents fixed and random terms separately rather than
only through dense observation covariance matrices:

```text
y = X beta + Z b + e
b ~ N(0, Lambda(theta) Lambda(theta)')
```

Grouped `Z`, independent-term variance estimation, Henderson prediction, and
estimated correlated blocks are now sparse. Correlated blocks use relative
Cholesky factors and a profiled residual scale, with one retained symbolic
factorization across covariance evaluations. Remaining pieces include more
complete boundary inference and selected sparse conditional covariances. Sparse CPU
factorization is the primary path. Dense or sufficiently
large frontal work is routed through JDistlib to GPU, oneMKL, OpenBLAS, or
portable CPU according to the existing policy.

## Pedigree performance

Pedigree prediction and the sparse variance fitter incorporate `A^-1`
directly. Multiple named pedigree terms, ordinary terms in the same model, and
unphenotyped ancestors are supported. Remaining work is scalable animal-level
PEV/reliability extraction and formula-level pedigree mapping. The dense
reference path continues to provide complete PEV/reliability results.

## GLMM likelihoods

PQL remains an explicitly named fast approximation. `glmer`-class parity
requires marginal maximum likelihood with a Laplace approximation, followed by
adaptive Gauss-Hermite quadrature for supported low-dimensional integrations.
These estimators must not be labeled REML: their likelihood and inference are
distinct from the current working-Gaussian REML calculation.

Accordingly, JLinAlg does not claim full `lme4` compatibility. Sparse Gaussian
unstructured covariance estimation, model-derived Satterthwaite/Kenward-Roger
inference, and formula-native constrained ML profiles now have R reference
tests, as detailed below. GLMM likelihoods are a separate implementation and
validation area; these Gaussian results do not certify them. `pedigreemm`'s core Gaussian animal model,
pedigree PQL facade, sparse `A^-1`, multiple pedigree terms, BLUP, and dense
PEV/reliability paths are present; scalable sparse PEV diagonals and
formula-native pedigree/new-data mapping remain open.

## User-facing compatibility

The compiled formula layer now covers fixed effects, contrasts, interactions,
offsets/weights, numeric and categorical random slopes, `||`, nested grouping,
and dense/sparse correlated blocks. Categorical slopes support treatment/sum
contrasts, full indicators, and `:`/`*` interaction expansion. As in lme4,
`||` splits formula terms, not every factor indicator column: a factor's
indicator columns remain a correlated block. Remaining formula work includes automatic
complete-case alignment for missing grouping rows and pedigree mappings.
Post-fit work still includes richer singular-fit diagnostics and boundary
profile coverage. Conditional/marginal prediction, response refit,
response simulation, and parametric bootstrap are now available as matrix-
first APIs; formula-native `newdata` compilation remains to be added.

## September 2026 Gaussian mixed-model repair and validation

`SparseUnstructuredCorrelatedModel.fit` estimates every lower-Cholesky entry
of every grouped block. The equation design is `Z L`, and the latent random
coefficients have covariance `sigma^2 I`. Only `sigma^2` is profiled; there is
no additional freely estimated block scale. Zero non-leading Cholesky
diagonals allow singular covariance fits while preserving positive marginal
variance bounds. Optimization uses bounded BOBYQA followed by a deterministic
pattern refinement and a local stationarity check. A singular result is
checked against two distinct interior starts: coordinate-wise stationarity
at a zero Cholesky diagonal can otherwise hide a better full-rank solution.
This repairs the categorical sum-contrast counterexample below, but is not a
guarantee of global optimality. The independent sparse
fitter also checks local stationarity instead of inferring success from an
evaluation count. The dense correlated reference optimizer now starts with
BOBYQA before pattern refinement.

`varianceBounds(minimum, maximum)` now constrains **physical** residual and
marginal random variances, not just relative-factor search coordinates. For
each relative marginal variance `r`, the residual-scale profile intersects
`[minimum/r, maximum/r]` with the residual variance bounds; the profiled scale
is clipped to that feasible interval. For unstructured blocks, `r` is the
corresponding squared Cholesky-row norm, not the squared diagonal alone.
Infeasible covariance candidates are rejected. These are variance bounds,
not independent box constraints on off-diagonal covariance entries.

The sparse likelihood centers the response within the fixed-effect mean
space and evaluates its quadratic as the conditional residual sum of squares
plus the sparse coefficient-precision penalty. It no longer subtracts large
`y'y`-sized quantities. This preserves fixed nonidentity/pedigree precision
semantics and prevents response-shift cancellation without allocating an
observation-square matrix.

`MixedFormula.compile(formula, table, FormulaOptions)` accepts residual
precision weights. All fixed and random designs and the offset-adjusted
response are whitened by `sqrt(weight)`; the response likelihood includes
the `sum(log(weight))/2` Jacobian. Fitted values, residuals, and offsets are
restored to observation units. Tests cover a correlated block plus a crossed
intercept, two separate correlated blocks, and correlated slopes without a
random intercept. The parser no longer consumes `offset(...)` as part of a
following random term.

Use `fitSparseUnstructured(...).correlatedFit()` or `fitCorrelated(...)` for
original-coordinate random modes and full covariance matrices. The wrapper's
`fit()` exposes sparse equation diagnostics and explicitly named **latent**
coefficients, not the original `b = L u` coefficients. `covarianceShape()` is
the first block's covariance relative to residual variance; full absolute
covariances are in `randomEffects()`.

Finite-DF inference is computed from the fitted model rather than supplied
derivatives. Central differences use **absolute, linear covariance-entry
coordinates**, including residual variance. Satterthwaite uses inverse
observed REML information and derivatives of `C = Var(beta)`. Kenward-Roger
uses expected information `-1/2 d2 log|V|_REML` and
`C_adjusted = C - sum(W_ij * d2C/dtheta_i dtheta_j)`. In linear covariance
coordinates this is the usual KR `P/Q` covariance adjustment. For the
coefficient-level, one-DF tests, KR moment matching reduces to `2/A2` using
the unadjusted contrast variance and expected information. This agrees with
the [pbkrtest covariance adjustment](https://github.com/hojsgaard/pbkrtest/blob/master/R/KR_vcovAdj.R)
and [coefficient-test moment calculation](https://github.com/hojsgaard/pbkrtest/blob/master/R/KR_modcomp.R).
No observation covariance or projection matrix is constructed. The original
sparse precision path, including nonidentity coefficient precision, remains
intact and is tested against an equivalent dense covariance model.

Formula-native `profileFixedEffect`, `profileResidualSd`, `profileRandomSd`,
and `profileCorrelation` re-optimize nuisance parameters under ML, even if
the supplied options request REML. Fixed coefficients are eliminated
algebraically, including the one-fixed-coefficient case. Random SDs are
constrained by a fixed Cholesky-row norm; two-coefficient correlations use
a fixed row angle. Nonconverged/nonfinite refits raise an error rather than
producing an apparently valid crossing.
Correlation profiles centered at `rho = -1` or `rho = 1` now evaluate the
estimable side and retain the boundary as the opposite endpoint. A boundary
endpoint is not an observed likelihood-ratio crossing: its corresponding
`lowerFound`/`upperFound` flag is false. The cutoff remains the regular
one-parameter chi-square cutoff, not a boundary-coverage correction.

### Frozen R accuracy evidence

`src/test/resources/r-reference/generate-mixed-sparse-reference.R` generated
four fixtures using R 4.6.1, lme4 2.0.6, lmerTest 3.2.1, and pbkrtest 0.5.5:
sleepstudy (180 rows), unbalanced sleepstudy (148), weighted crossed effects
(240), and two unstructured blocks (240). The weighted KR reference uses
explicitly whitened designs, preserving the same likelihood/estimand while
giving pbkrtest its identity residual basis.

Across these fixtures, log-likelihood discrepancies are below `2e-11`,
maximum absolute covariance-entry discrepancy is below `0.00085`, and
coefficient SE discrepancies are below `4e-6`. Finite DF differences are
below `0.0006`. The unbalanced model has a nonzero KR adjustment: its
intercept SE changes from approximately `6.81795` to `6.86519`, matching R.
The independent balanced-model test recovers the ANOVA denominator DF, and
the dense/sparse test compares covariance estimates, BLUPs, and fitted values.

| Sleepstudy ML 95% profile | Java lower | Java upper | R lower | R upper |
|---|---:|---:|---:|---:|
| Fixed slope | 7.358652 | 13.575920 | 7.358653 | 13.575919 |
| Random intercept SD | 14.381384 | 37.716055 | 14.381474 | 37.715996 |
| Random slope SD | 3.801177 | 8.753392 | 3.801164 | 8.753383 |
| Residual SD | 22.898256 | 28.857999 | 22.898267 | 28.857997 |
| Random correlation | -0.481493 | 0.685019 | -0.481501 | 0.684986 |

All 34 isolated mixed/formula tests passed after the counterexample and
categorical-slope repairs. Source and tests were compiled
with `javac` against the existing build output in `build/mixed-classes` and
`build/mixed-test-classes`; the coordinated parent build owns whole-repository
Gradle verification.

### Reviewer counterexamples and categorical formula evidence

`src/test/resources/r-reference/generate-mixed-review-reference.R` freezes
the three reviewer reproductions. `SparseMixedReviewRegressionTest` checks:

- **Large response shift:** for the 40-row grouped-intercept example, adding
  `1e8` changes fitted variances by less than `4e-8`. The shifted independent
  fit returns random variance `1.944899183164`, residual variance
  `0.781028127606`, and REML log-likelihood `-61.460390946826`; R's unshifted
  values are `1.944899152242`, `0.781028128338`, and `-61.460390910556`.
  Both independent and unstructured fitters pass the shift regression.
- **Physical bounds:** `[.01, .1]` returns both variances at `.1` and REML
  log-likelihood `-185.300514613339` in both fitters. The R reference
  independently optimizes the same absolute-variance box likelihood, rather
  than treating lme4's relative Cholesky bounds as physical bounds. A separate
  correlated-slope test checks both marginal variances and residual variance.
  Finite-DF requests at variance bounds are rejected rather than assigned
  unstable, unconstrained-information DF.
- **Boundary correlation profiles:** the positive-correlation example gives
  `[0.996388966839, 1]`; reversing the random-slope sign gives
  `[-1, -0.996388966839]`. Independently constrained R/lme4 profile refits give
  the interior endpoint magnitude `0.996388964706`, a difference of about
  `2.1e-9`. The available-side crossing flags are true and boundary-side flags
  are false.

`src/test/resources/r-reference/generate-mixed-categorical-reference.R`
generates **40 random-design comparisons** (20 expressions under each of
treatment and sum coding) and **five fitted-model comparisons** against
lme4. `MixedCategoricalFormulaTest` checks column names, values, interaction
ordering, marginal-term contrast selection, no-intercept indicators, and
`||` term splitting. Numeric controls remain covered by the existing tests.
Levels use the table's encounter order; R factors are explicitly assigned the
same level order. Missing categorical slopes are rejected, not silently
omitted.

The five 600-row, 50-group fits use the same precision weights and offsets in
Java and R, and compare covariance estimates, fixed coefficients, random
modes, and observation-unit fitted values. Log-likelihood discrepancies are
below `3e-11`. Regression gates use `1e-3 * max(1, abs(reference))` for
covariance entries and `5e-4 * max(1, abs(reference))` for modes/fitted values;
these are test tolerances, not claims of measured maximum error. The sum-coded
case verifies that the restart check escapes the poorer singular solution.

| Random expression / coding | Java log-likelihood | Java fit seconds | R lme4 fit seconds |
|---|---:|---:|---:|
| `1+f`, treatment | -981.680378339101 | 0.112 | 0.050 |
| `0+f`, full indicators | -981.680378339101 | 0.070 | 0.010 |
| `1+f`, sum | -981.680378339102 | 0.189 | 0.030 |
| `0+f:x`, interaction | -1159.028797782225 | 2.638 | 0.040 |
| `1+f*x`, factorial | -725.222757842871 | 0.877 | 0.190 |

These categorical timings are single-run fit measurements, not benchmark
medians; they exclude compilation and R package loading and use
residual-approximation DF in Java. The interaction fit is genuinely singular
in both implementations and is especially slow in Java (20,699 objective
evaluations including restarts). Java is slower than R on all five measured
categorical cases. They do not establish categorical finite-DF parity at a
singular estimate.

### Measured performance and memory

`SparseMixedParityBenchmark` consumes fit results and gates timing output on
convergence, likelihood accuracy, and reference DF. Java figures below are
median seconds of three fits after one warm-up, portable CPU, `-Xmx96m`.
R fitting figures include lmerTest covariance-derivative construction;
the separate R KR column is inference on an already fitted equivalent model.
Consequently Java KR includes fitting work that the R KR column excludes.
Small R times have coarse timer resolution and all results are host-specific.

| Fixture | Java Satterthwaite fit | Java KR fit | R lmerTest fit | R KR inference only |
|---|---:|---:|---:|---:|
| Sleepstudy, 180 rows | 0.017861 | 0.010438 | 0.010 | 0.100 |
| Unbalanced, 148 rows | 0.005858 | 0.007273 | 0.030 | 0.050 |
| Weighted crossed, 240 rows | 0.058977 | 0.052707 | 0.030 | 0.330 |
| Balanced, 12,000 rows | 0.846031 | 0.845339 | 0.220 | Not run |

Java measurements above were refreshed after the stable-quadratic,
physical-bound, and singular-restart solver changes. R columns retain their
frozen reference-run timings; they are not simultaneous paired timings.

`src/benchmark/r/mixed_sparse_scalability.R` generates the 12,000-row fixture
and R timing. Java holds 3,000 random coefficients with 4,500 equation and
factor nonzeros and completes both inference methods in a **96 MB heap**.
A single dense observation matrix would require 1.15 GB. Large-model
Satterthwaite DF are `1499.002799` and `1499.003827` versus R's 1499; KR DF are
`1499.000202` and `1498.999425`, within `0.0006` of the balanced-model identity
1499. The finite-difference
approximation is visible at this scale. R KR was not run for this large
fixture because pbkrtest constructs dense observation-scale matrices.

Java is not uniformly faster than R: crossed-block and large-model
Satterthwaite remain slower. Reusing symbolic factorization reduced the
initial 12,000-row Java timing from roughly 13 seconds to under one second.

### Explicit remaining limits

- These are Gaussian ML/REML results, not a full lme4/pedigreemm or GLMM
  compatibility claim.
- Finite DF require REML and a converged, identifiable **interior** covariance
  estimate. Boundary fits can be returned with residual-approximation DF;
  singular finite-DF information and estimates numerically near physical
  variance bounds are rejected. Joint multi-DF KR F tests and
  boundary-adjusted finite-DF inference are not implemented by this API.
- The derivative machinery uses scaled central differences, not analytic
  covariance derivatives. It takes quadratically many sparse evaluations in
  covariance-parameter count and stores parameter/fixed-effect-sized arrays.
  Sparse factor fill can still grow substantially for crossed designs.
- Profiles currently use regular one-parameter chi-square cutoffs. Automatic
  profiles centered at zero random SD, general correlations in blocks larger
  than two coefficients, and nonregular boundary-coverage corrections remain
  open. Automatic handling of SD profiles reaching physical variance-box
  boundaries is also incomplete; infeasible constrained refits raise errors.
  The supported correlation profiles centered at `+/-1` retain that boundary
  with its crossing flag false; unfound general bounds remain NaN.
- Full original-coordinate random-effect conditional covariance/PEV blocks
  are not exported by the unstructured result. Its latent equation PEVs must
  not be interpreted as original-coordinate PEVs.
- The new grouped unstructured optimizer does not yet jointly estimate a
  grouped covariance and arbitrary pedigree precision scales in one model.
  The existing multiple-precision sparse/pedigree path remains supported;
  formula-level pedigree mapping and new-data compilation remain open.
- The specialized changing-term scan fitter retains residual-approximation
  DF. Ordinary prepared sparse fits now support Satterthwaite and KR.
- Backend calls support the existing prepared native sparse-factor interface,
  but the measured validation here used portable CPU. No CHOLMOD native run
  was certified on this host. Sparse `PREFERRED` uses the existing sparse
  backend preference, avoiding repeated GPU setup for host sparse equations.
- Automatic complete-case alignment, arbitrary R formula expressions, and
  complete reproduction of lme4 optimization controls remain incomplete.
  Categorical treatment/sum slopes, full indicators, and `:`/`*` interactions
  are implemented and tested as described above. Physical marginal variance
  bounds are enforced, but neither the local stationarity check nor the
  singular-result restarts guarantee a global covariance optimum.
