# GAM, GAMM, and distributional-regression roadmap

This roadmap makes one basis-and-penalty representation the boundary between
model construction and numerical fitting. A smooth contributes a row-major
basis `B`, one or more positive-semidefinite penalties `S_j`, identifiability
constraints, prediction metadata, and optional sparsity metadata. GAM, GAMM,
and distributional models therefore share basis construction and differ only
in their likelihood and covariance layers.

The implementation order prioritizes reusable compilation and repeated-fit
speed for GWAS, TWAS, EWAS, and PWAS. GPU execution uses JDistlib automatic
routing when available; the CPU fallback remains oneMKL, OpenBLAS, and then the
portable backend.

## Stage 1: penalized GAM core

### Model representation

- Treat every smooth as basis plus penalty, never as a callback executed during
  each solver iteration.
- Center ordinary smooths at compilation time and move every penalty null space
  into the unpenalized fixed-effect design.
- Allow multiple penalties per basis from the beginning. This is required for
  tensor products, adaptive smooths, and later GAMMs.
- Retain knots, transforms, constraints, and factor levels so prediction never
  rebuilds training-time choices.

### Initial bases

1. P-splines with first- and second-difference penalties.
2. Cubic regression splines and cyclic cubic splines.
3. Thin-plate regression splines for `mgcv` compatibility.
4. Tensor products `te()` and interaction-only `ti()` constructed from marginal
   bases and Kronecker penalties.
5. Factor-by, varying-coefficient, random-effect, and Markov-random-field
   smooths.

### Fitting and smoothing selection

- Exact Gaussian REML is the first reference path. A single-penalty smooth is
  decomposed into an unpenalized null space and a whitened random-effect design;
  the smoothing parameter is residual variance divided by smooth variance.
- Penalized IRLS adds binomial, Poisson, gamma, inverse-Gaussian, negative
  binomial, Tweedie, and other existing GLM families.
- REML is the default. ML, GCV/UBRE, fixed smoothing parameters, and shared or
  linked smoothing parameters are explicit options.
- Results report parametric beta/SE/statistic/p-value, term EDF, smoothing
  parameters, term contributions, omnibus smooth tests, covariance, intervals,
  deviance diagnostics, and prediction.

### Large-data execution

- `PreparedGam` compiles marginal bases, constraints, and penalty factorizations
  once for repeated responses.
- A BAM-like mode forms cross-products in blocks without materializing the full
  design. Optional covariate discretization stores unique marginal basis rows
  plus integer row maps.
- Batched dense products and factorizations route through JDistlib. Host-side
  model compilation remains parallel, but nested BLAS/GPU oversubscription is
  avoided.
- Association scans residualize the phenotype, tested marker, or both against a
  prepared smooth null model. The nuisance smooth is not refitted for every
  marker unless the requested test requires it.

## Stage 2: unified GAMM

The predictor is `eta = X beta + sum_j B_j alpha_j + Z b`, where each smooth
penalty and random-effect precision is a named block. The covariance builder
accepts ordinary grouped effects, correlated slopes, pedigree relationship
matrices, genomic relationship matrices, cryptic relatedness, sparse precision
matrices, and AR/ARMA/ARIMA residual correlation.

- Gaussian responses use exact ML/REML with dense and sparse paths.
- Non-Gaussian responses retain PQL as a fast compatibility path and add a
  Laplace marginal-likelihood path as the default accuracy path.
- Multiple penalties attached to one smooth remain first-class; the design must
  not inherit `gamm4`'s `te()` and adaptive-smooth restrictions.
- Residual-approximation, Satterthwaite, and Kenward-Roger remain available for
  parametric terms where they apply.
- Results add simultaneous and pointwise smooth intervals, conditional and
  marginal prediction, random/smooth modes, PEVs, and new-level policies.

## Stage 3: distributional and vector additive models

A family declares one linked additive predictor for each modeled parameter,
not hard-coded `mu`, `sigma`, `nu`, and `tau` fields. This covers GAMLSS-style
location/scale/shape models and VGAM-style vector predictors with one solver
contract. Predictors may have separate or shared smooths and linear equality or
proportionality constraints.

Initial families are Gaussian location-scale, gamma mean/shape, beta
mean/precision, negative-binomial mean/dispersion, zero-inflated and hurdle
counts, ordinal and multinomial responses, followed by selected flexible
four-parameter continuous families. Each family supplies parameter links, valid
domains, log likelihood, score, and expected or observed information.

Block Fisher scoring or Newton updates all distribution parameters, with inner
penalized weighted solves and outer smoothing-parameter optimization. Step
halving, trust-region bounds, domain transforms, and warm starts are mandatory.
Results include quantile residuals, centile prediction, parameter-specific term
summaries, and likelihood-based model comparison.

## Validation gates

Every milestone requires deterministic fixtures generated in R, using `mgcv`
first, `gamlss` for distributional models, `gamm4`/`lme4` for mixed models, and
VGAM for vector responses. Tests compare basis values, penalties, smoothing
parameters, EDF, coefficients, fitted values, likelihood/deviance, predictions,
and term tests. Performance claims require separate cold-start, prepared-model,
large-`n`, and association-scan benchmarks.

## Completion order

1. Deliver the single-penalty Gaussian correctness path and R fixtures.
2. Add penalized IRLS and fixed/GCV/ML smoothing selection.
3. Add multi-penalty and tensor-product bases.
4. Add `PreparedGam`, block/discrete basis execution, and association scans.
5. Compose smooth blocks with existing LMM, pedigree, GRM, and residual models.
6. Add Laplace GAMM and validate against `gamm4`/`lme4`.
7. Introduce the multi-predictor family interface and first distributional
   families.
8. Add VGAM-style multinomial/ordinal constraints and vector responses.
