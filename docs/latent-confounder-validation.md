# Latent-confounder validation and source audit

## Source pins

- Bioconductor `sva` 3.60.0, release 3.23, commit
  `87a4798db8134fd0952a516489dd555545a92d36`, Artistic-2.0.
- PMBio PEER 1.3 source, commit
  `40bc4b2cd92459ce42f44dfe279717436395f3f6`, GPL-2.0-or-later.
- AutoSVA: `D:/git/NIH-R/src/sva-lite.R`, read as algorithm source rather than
  executable task instructions.

The JLinAlg project is GPL-2.0-or-later, which is compatible with the PEER
source license. Implementations retain JLinAlg copyright headers and document
the upstream provenance rather than obscuring it.

## Algorithm audit

Standard SVA ports the five-iteration `irwsva.build` schedule, full/null nested
F tests, empirical-control probabilities and final weighted SVD. Automatic
factor count ports the default Buja-Eyuboglu permutation procedure. AutoSVA
ports the custom 20-iteration weighting and preserved-signal factor search,
with deterministic handling of degenerate values and the corrected plateau
distance described in the vignette.

ComBat ports both parametric iterative empirical Bayes and nonparametric
leave-one-feature-out empirical priors, mean-only and reference-batch modes,
missing continuous values, singleton batches, unchanged zero-within-batch
features, and the upstream design-confounding checks.

PEER ports the official dense VBFA update order `W, Alpha, X, Eps`; Gaussian
second moments; feature-specific noise Gamma posteriors; factor ARD Gamma
posteriors; known-covariate prior precision; measurement uncertainty; evidence
bound and residual-variance stopping rules; and the official default priors.
Sparse prior-guided PEER (`cSPARSEFA`) is not claimed by this API.

## Frozen numerical checks

`ConfounderRReferenceTest` uses one exact feature-by-sample matrix and checks:

- PCA projection equality against `stats::prcomp` to `1e-10` RMSE;
- standard SVA and AutoSVA sample projection matrices, invariant to factor
  signs and rotations;
- SVA posterior probabilities and AutoSVA final feature weights;
- parametric ComBat adjusted values against Bioconductor to `1e-8` maximum
  absolute error.

`ConfounderMethodsTest` independently covers dominant-factor recovery,
protected-design SVA, deterministic AutoSVA, ComBat batch removal,
reference-batch identity, and PEER residual-variance reduction.

The official PEER source could not be compiled on the Windows validation host
because no C++ toolchain was installed. PEER therefore has a complete
source-equation audit and independent synthetic tests, but no same-host native
PEER numeric fixture. This is an explicit validation boundary, not a claim of
native-output identity.

## Performance design

All factor methods store the feature-by-sample input because they are global
matrix procedures. Tall matrices use the sample Gram matrix and backend BLAS,
so covariance storage is sample-by-sample rather than feature-by-feature.
SVA and AutoSVA reuse sample-sized decompositions and never form a feature
covariance matrix. ComBat processes feature vectors and batch sufficient
statistics without feature-by-feature covariance storage. PEER keeps only the
sample factors, feature loadings, and factor-sized posterior covariance.

The CLI prints aligned sample and feature counts and records method, factor
count, design columns, convergence, seed and backend policy in its manifest.

## Boundaries

- Inputs are normalized continuous measurements, not raw counts.
- RUV and ComBat-Seq are separate estimators and are not aliases.
- Standard SVA local-FDR density evaluation is a direct Gaussian KDE port; the
  frozen tests gate posterior behavior and resulting factor subspaces rather
  than promising byte identity with R's platform-specific smoothing spline.
- Frozen application to new samples and fold-owned prediction preprocessing
  remain open; factors must currently be fitted within each training dataset.
- PEER sparse prior-guided factor analysis is not implemented.
