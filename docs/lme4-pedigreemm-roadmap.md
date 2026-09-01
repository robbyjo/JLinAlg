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
- Sparse pedigree variance estimation using `A^-1` directly.
- Cholesky-parameterized unstructured correlated random blocks on the dense
  reference likelihood path.
- Formula `||` independent terms, nested grouping shorthand, correlated
  single-bar blocks, and `fixef`/`ranef`/`VarCorr`-style accessors.
- ML-only nested likelihood-ratio model comparison.
- Retained-design response refits with variance-component warm starts.
- Marginal and conditional prediction on new data, including an explicit
  zero-mode policy for unseen grouping levels.
- Multiple named pedigree precision terms and ordinary independent random
  terms in the same sparse REML model, including unphenotyped ancestors.

## High-performance mixed-model core

The next engine must represent fixed and random terms separately rather than
only through dense observation covariance matrices:

```text
y = X beta + Z b + e
b ~ N(0, Lambda(theta) Lambda(theta)')
```

Grouped `Z`, independent-term variance estimation, and Henderson prediction
are now sparse. Remaining pieces are moving correlated blocks onto the sparse
precision path, full boundary-aware optimization, and selected sparse
conditional covariances. Sparse CPU
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

Accordingly, JLinAlg does not yet claim full `lme4` likelihood parity. The
largest remaining parity items are sparse correlated-block likelihoods,
Laplace/AGQ GLMM estimation, boundary-aware optimizer behavior, and
profile/bootstrap inference. `pedigreemm`'s core Gaussian animal model,
pedigree PQL facade, sparse `A^-1`, multiple pedigree terms, BLUP, and dense
PEV/reliability paths are present; scalable sparse PEV diagonals and
formula-native pedigree/new-data mapping remain open.

## User-facing compatibility

The compiled formula layer now covers fixed effects, contrasts, interactions,
offsets/weights, independent random intercepts/slopes, `||`, nested grouping,
and dense correlated blocks. Remaining formula work includes automatic
complete-case alignment for missing grouping rows and pedigree mappings.
Post-fit work still includes response simulation, richer singular-fit
diagnostics, profile likelihood, and parametric bootstrap. Conditional and
marginal prediction plus response refit are now available as matrix-first APIs;
formula-native `newdata` compilation remains to be added.
