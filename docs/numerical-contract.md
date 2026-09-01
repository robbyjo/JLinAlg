# Numerical contract

## Storage and validation

Internal dense matrices are contiguous row-major `double[]` arrays, matching
JDistlib's unified FP64 algebra contract. Public two-dimensional-array overloads
validate rectangular shape and copy into that representation. Responses,
designs, covariance bases, and starting values must be finite. Covariance bases
must be symmetric to a scale-aware tolerance.

## OLS

For a full-column-rank design, JLinAlg solves

```text
minimize ||y - X beta||²
```

using column-pivoted QR. It does not form the normal equations to estimate
`beta`. The coefficient covariance is reconstructed from the triangular QR
factor and pivot, then scaled by

```text
s² = residual sum of squares / (observations - rank).
```

Inference uses a Student t distribution with `observations - rank` degrees of
freedom. The reported likelihood is the maximized ordinary Gaussian likelihood,
which uses `RSS / observations` as its variance estimate.

When explicitly enabled, the rank-deficient path uses an SVD tolerance of

```text
max(rows, columns) * ulp(1.0) * largest singular value.
```

It returns the Moore-Penrose minimum-norm coefficient vector and generalized
covariance. This is a numerical convention, not proof that individual original
coefficients are estimable.

## Gaussian penalized regression

For weights normalized to sum to `n`, JLinAlg minimizes

```text
(1 / 2n) sum_i weight[i] (y[i] - intercept - x[i]' beta)^2
+ lambda ((1 - alpha) / 2 sum_j factor[j] beta[j]^2
          + alpha sum_j factor[j] |beta[j]|).
```

`alpha = 0` is ridge, `alpha = 1` is LASSO, and intermediate values are elastic
net. Cyclic coordinate descent uses residual updates and warm starts along a
strictly descending lambda path. Unless disabled, weighted means and weighted
root-mean-square scales standardize each predictor for fitting; results are
transformed back to the caller's original predictor scale. A constant predictor
has scale one and remains uninformative. The intercept is never penalized.

The automatic lambda maximum is the first L1 threshold at which all penalized
slopes are zero. Pure ridge has no finite all-zero threshold, so automatic ridge
paths use a small effective alpha only to choose a useful scale; every fitted
model still uses `alpha = 0`. Cross-validation learns centering and scaling from
each training fold and reports weighted validation mean squared error.

Ridge inference uses the effective model dimension

```text
df_model = intercept_dimension + trace((G + lambda D)^-1 G)
```

and `df_residual = n - df_model`, where `G = X' W X / n` after preprocessing
and `D` contains penalty factors. Its covariance is the model-based Gaussian
sandwich covariance for the ridge estimator. LASSO/elastic-net active-set OLS
refits are explicitly conditional on selection; they are not selective
inference and do not correct for the search over predictors.

## REML

For `V = sum variance[k] K[k]`, the restricted log likelihood is

```text
-0.5 * ((n-p) log(2 pi)
        + log|V|
        + log|X' V^-1 X|
        + y' P y)
```

where

```text
P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1.
```

Both `V` and `X' V^-1 X` must be positive definite. Fixed effects therefore
must have full column rank under the fitted covariance. Variances are optimized
as logarithms, preserving positivity. The analytical log-scale score and
expected Fisher information are

```text
score[k] = 0.5 variance[k]
           * (y' P K[k] P y - trace(P K[k]))

information[k,l] = 0.5 variance[k] variance[l]
                   * trace(P K[k] P K[l]).
```

Fisher steps are regularized, limited, and backtracked until the restricted
likelihood increases. A normalized gradient step is tried if the Fisher step
cannot improve it. Convergence requires the projected score tolerance; a
likelihood-change tolerance is additionally used after an accepted step.
Variance bounds are explicit options, and the projected score accounts for an
active bound.

The present implementation materializes dense `n x n` matrices, including
`P`, so its memory complexity is `O(n²)` and its leading runtime is `O(n³)`.
It is suitable as a clear reference implementation and for moderate sample
sizes.

### Fixed-effect association inference

For every fixed-effect coefficient, REML reports beta, its model-based standard
error, a Wald t statistic, denominator degrees of freedom, and a two-sided
p-value. The default denominator DF is the constant fast approximation

```text
df = n - rank(X) - 1.
```

With `DegreesOfFreedomMethod.SATTERTHWAITE`, let `C = vcov(beta)`, let `g_j` be
the gradient of `C[j,j]` with respect to the log variance parameters, and let
`S` be their covariance matrix obtained by inverting the expected REML
information. JLinAlg uses

```text
df[j] = 2 C[j,j]^2 / (g_j' S g_j).
```

The covariance derivatives are analytic for the covariance-component model.
The expected Fisher information is already produced by fitting, so the option
does not require coefficient-by-coefficient refits. A small numerical ridge is
used only when factoring the variance-parameter information matrix.

With `DegreesOfFreedomMethod.KENWARD_ROGER`, JLinAlg implements the original
Kenward-Roger calculation for covariance structures that are linear in the raw
variance components. Let `Phi = (X' V^-1 X)^-1`, let `W` be the covariance of
the raw REML variance estimates, and define

```text
P[i]    = -X' V^-1 K[i] V^-1 X
Q[i,j]  =  X' V^-1 K[i] V^-1 K[j] V^-1 X
U       = sum_i sum_j W[i,j] (Q[i,j] - P[i] Phi P[j])
Phi_A   = Phi + 2 Phi U Phi.
```

`Phi_A` supplies the reported fixed-effect SE. For each one-dimensional
coefficient contrast, denominator DF follows the Kenward-Roger moment-matching
calculation using `P`, `W`, and the unadjusted `Phi`. The corresponding
one-dimensional F scaling is exactly one, so the reported statistic is the
usual signed t statistic based on `Phi_A` with that denominator DF.

The fitted covariance bases are linear in their raw nonnegative variance
components, so second derivatives of `V` with respect to those parameters are
zero. Known covariance contributions are included in `V` but correctly omitted
from `W` because they were not estimated. The implementation follows the dense
reference route and is therefore intentionally an opt-in, relatively expensive
calculation.

### Gaussian mixed-model simulation and bootstrap

Marginal LMM simulation uses the fitted fixed mean and independently draws
each scalar random term and residual:

```text
y* = X beta_hat + sum_j Z_j b_j* + e*
b_j* ~ N(0, sigma_j^2 I)
e*   ~ N(0, sigma_e^2 I).
```

Conditional simulation replaces each `b_j*` with its fitted conditional mode
and redraws only `e*`. Pedigree marginal simulation instead draws the complete
breeding-value vector from `N(0, sigma_a^2 A)` by Cholesky factorization of the
numerator relationship matrix; conditional simulation retains the fitted
BLUPs.

Parametric bootstrap always uses marginal simulation, then warm-refits the
same retained design. Replicate seeds are derived independently from the
requested seed and replicate index, so results do not depend on CPU scheduling
order. Summaries use the successful replicates to report empirical bias,
sample standard deviation, and linearly interpolated percentile intervals.
Nonconverged and failed refits are excluded from those summaries and retained
as structured failures; they are never silently replaced.

## Pedigree REML

The pedigree layer constructs the numerator relationship matrix `A` by the
tabular recurrence. Input order need not place parents before offspring. Every
named parent must have its own pedigree entry; unknown parents are represented
by `null` and contribute zero relationship. Duplicate identifiers,
self-parenting, identical sire/dam identifiers, missing parent records, and
ancestry cycles are rejected.

Observation identifiers construct an incidence matrix `Z` without explicitly
materializing it. Pedigree REML calls the dense Gaussian engine with covariance

```text
V = additiveVariance Z A Z' + residualVariance I.
```

Breeding values use `additiveVariance A Z' V^-1 (y - X beta)`. Prediction-error
variances are the diagonal of

```text
additiveVariance A - G Z' P Z G,
```

with variance components treated as estimated constants. Reliabilities divide
that prediction-error variance by each individual's prior additive variance.
The original implementation is a dense reference path: memory grows as
`O(max(n²,m²))` for `n` observations and `m` pedigree individuals. Large
pedigrees will require sparse `A^-1` and Henderson mixed-model equations.

`Pedigree.sparseRelationshipMatrixInverse()` already constructs `A^-1`
directly as full symmetric CSR using `A = T D T'` and
`A^-1 = T'^-1 D^-1 T^-1`. Each individual's Mendelian-sampling row has at most
three nonzero coefficients (self, sire, dam), so construction avoids a dense
matrix inversion. `SparsePedigreeReml` integrates this CSR precision into
sparse REML/ML variance estimation and BLUP equations without materializing
`A` or `Z A Z'`. Its current high-throughput result omits full animal-level
PEV/reliability diagonals; use the dense reference fitter when those quantities
are required.

## ARIMA and seasonal ARIMA

After applying ordinary and seasonal differences, the stationary series uses
the R-compatible equation

```text
z[t] - mu = sum_i phi[i] (z[t-i] - mu)
            + epsilon[t] + sum_j theta[j] epsilon[t-j].
```

Nonseasonal and seasonal AR and MA polynomials are multiplied, including their
interaction lags. Optimization variables map through partial autocorrelations;
AR coefficients are stationary and MA coefficients are invertible. Conditional
innovations before the largest effective AR lag initialize past innovations to
zero and are omitted from the objective. For `nUsed` retained innovations,

```text
sigma2 = sum epsilon[t]^2 / nUsed
logLik = -nUsed/2 * (log(2 pi) + 1 + log(sigma2)).
```

The reported likelihood and information criteria are conditional, not the exact
state-space likelihood reported by an exact-ML ARIMA fit. Forecast means recurse
through the fitted ARMA and inverse-differencing polynomials. Forecast variances
use the MA-infinity impulse weights, convolved with the inverse-differencing
filter for integrated models.

Stationary ARMA correlation matrices use

```text
gamma[k] = sum_j psi[j] psi[j+k]
correlation[i,j] = gamma[abs(i-j)] / gamma[0],
```

with adaptive impulse truncation. An AR(1) therefore yields the exact Toeplitz
sequence `rho^abs(i-j)` to numerical precision.

## LMMs with ARIMA errors

For an undifferenced error order, the mixed-model covariance is

```text
V = sum_k variance[k] Z[k] Z[k]' + residualVariance R_ARMA(theta).
```

ARMA parameters are transformed to the stationary/invertible region. For each
candidate, the variance components and fixed effects are fitted by REML; the
outer derivative-free optimization profiles that restricted likelihood over
the ARMA parameters. This is computationally more expensive than an ordinary
LMM, but the final fixed-effect inference still supports the fast residual-DF,
Satterthwaite, and Kenward-Roger policies.

For nonzero ordinary or seasonal integration, the differencing polynomial is
applied identically to `y`, `X`, and every `Z` before profile REML. This is a
regression/mixed model in differences with stationary ARMA errors. Columns
annihilated by differencing are rejected rather than silently dropped.

## GLM

For mean `mu`, link `g`, variance function `v`, prior weight `a`, and linear
predictor `eta = X beta + offset`, IRLS uses

```text
working response = eta + (y - mu) / (d mu / d eta) - offset
working weight   = a (d mu / d eta)^2 / v(mu).
```

The weighted design and response are solved by pivoted QR. Candidate
coefficient steps are halved until the model deviance does not increase.
Convergence requires both relative coefficient and relative deviance changes to
meet the configured tolerance.

The built-in contracts are:

- Gaussian variance `1` with identity link;
- binomial variance `mu(1-mu)` with logit link;
- Poisson variance `mu` with log link.

Gaussian dispersion is the Pearson estimate by default. Binomial and Poisson
use dispersion one unless the caller explicitly requests Pearson dispersion.
Inference is asymptotic normal/Wald inference. Binomial log likelihood and AIC
are `NaN` for nonintegral effective trial/success counts, while deviance and
quasi-likelihood estimation remain available.

## GLMM PQL

For fixed-dispersion non-Gaussian families, the PQL working residual covariance
is diagonal with entries

```text
R[i,i] = v(mu[i]) / (priorWeight[i] (d mu[i] / d eta[i])^2).
```

At every outer iteration, JLinAlg fits

```text
working response - offset ~ N(X beta, G + R)
G = sum variance[k] K[k]
```

by REML, treating `R` as known and estimating only the random covariance
scales. The conditional random linear predictor is the working-model BLUP

```text
G (G + R)^-1 (working response - offset - X beta).
```

Outer convergence requires relative changes in the total linear predictor and
all variance scales plus convergence of the final inner REML fit. The reported
working restricted likelihood is conditional on the final linearization; it is
not an exact marginal GLMM likelihood and must not be used as though likelihoods
from different approximations were directly comparable.

Fixed-effect beta, SE, t statistic, DF, and p-value are inherited from the
final working REML fit. Consequently, both the fast residual approximation and
the Satterthwaite option describe the final PQL linearization, not an exact
finite-sample distribution under the original non-Gaussian mixed model.

First-order PQL is known to be least reliable for binary outcomes with small
clusters, rare events, or substantial heterogeneity. These limitations are a
property of the estimator rather than the linear-algebra backend.

## Genomic relationship matrices and cryptic relatedness

For each retained variant `j`, called alternate-allele dosages estimate
`p_j = mean(g_j)/2`. Missing calls are mean-imputed, hence have centered value
zero. The standardized matrix and additive GRM are

```text
Z_ij = (g_ij - 2 p_j) / sqrt(2 p_j (1-p_j))
K = Z Z' / m,
```

where `m` is the number of variants passing the requested MAF and call-rate
filters. Matrix multiplication uses the selected JDistlib backend. The GRM is
intended to be built from suitably QC'd, preferably LD-pruned markers.
Off-diagonal additive relationship is twice the reported kinship coefficient.

`K` is a covariance basis for REML, GLMM PQL, P3D/EMMAX, and the retained-REML
Burden/SKAT/SKAT-O score null. Repeated observations use `Z_subject K
Z_subject'`. Cox kinship frailty uses `(K + epsilon mean(diag(K)) I)^-1` as
Gaussian precision, with caller-visible relative `epsilon`; this regularizes
duplicate samples and marker-rank deficiency. It remains a dense reference
path.

## Cox proportional hazards and Gaussian frailty

For event time `t`, the risk set follows counting-process convention

```text
R(t) = { i : start_i < t <= stop_i }.
```

The fixed model maximizes the Cox log partial likelihood. Breslow ties use one
risk-set denominator raised to the event count. Efron ties progressively remove
`k/d` of the tied-event risk sum for `k = 0,...,d-1`. Newton updates use the
analytic score and observed information; exponentials are shifted by the
largest linear predictor in each active risk set. Coefficient inference is an
asymptotic Wald z test, and hazard ratios plus confidence limits are obtained
by exponentiating beta and its normal-theory interval.

For right-censored input, risk moments are accumulated after one descending
time sort per stratum. General delayed-entry/start-stop input retains direct
risk-set evaluation. The reported baseline hazard uses the selected tie rule;
baseline survival is `exp(-cumulativeHazard)`.

For Gaussian random terms, conditional modes maximize

```text
l_partial(beta, b) - 1/2 sum_j b_j' P_j b_j / sigma_j^2.
```

Frailty variances are profiled on the log scale using the Laplace objective

```text
l_Laplace = l_penalized
            + 1/2 log|Q|
            - 1/2 log|H_bb|,
```

where `Q` is block-diagonal Gaussian precision and `H_bb` is the random-effect
block of penalized observed information at the conditional mode. Fixed-effect
covariance is the fixed block of the inverse joint penalized information.
Pedigree frailty uses the directly constructed numerator-relationship
precision `A^-1 / sigma_a^2`; unphenotyped ancestors remain in the coefficient
system. This is a dense Laplace reference likelihood for the random-effect
block, not REML, gamma frailty, or adaptive quadrature.

## Mendelian randomization

Summary associations are restricted initially to biallelic SNPs with alleles
`A`, `C`, `G`, or `T`, finite effects, and positive finite standard errors.
Harmonization aligns outcome effects to exposure effect alleles. Palindromic
variants are retained only when both effect-allele frequencies identify one
orientation, but not its reverse, within the configured frequency tolerance.

For independent variants, first-order IVW fits the zero-intercept weighted
regression

```text
betaY[j] = causalEffect betaX[j] + error[j]
weight[j] = 1 / seY[j]^2.
```

It treats exposure associations as fixed. The fixed-effect standard error is
`1 / sqrt(sum(weight betaX^2))`. Multiplicative random effects multiply its
variance by `max(1, Q/(J-1))`; underdispersion is not used to make inference
more precise. Cochran's Q and I-squared describe between-instrument
heterogeneity. Wald-ratio standard errors likewise use the first-order
`seY / abs(betaX)` approximation. Approximate F statistics are
`(betaX / seX)^2`.

MR-Egger is weighted regression with an intercept. Instruments are first
reoriented to positive exposure effects, flipping their outcome effects too.
Its standard errors use multiplicative overdispersion capped below at one.
The slope depends on the InSIDE assumption, and low I-squared GX warns that
regression dilution from imprecise exposure associations may be material.

The weighted median uses ratio weights `betaX^2 / seY^2` and linear
interpolation at cumulative probability 0.5. Its standard error is a
reproducible parametric bootstrap that independently samples exposure and
outcome associations from their reported normal sampling distributions. The
weighted median requires at least 50% of its weight to arise from valid
instruments for consistency.

For correlated variants, generalized IVW and MR-Egger use

```text
Omega = diag(seY) R diag(seY),
```

where `R` is a symmetric positive-definite LD correlation matrix aligned to the
reported effect alleles and instrument order. Cholesky solves are used rather
than an explicit inverse. Reorienting an MR-Egger variant changes the signs of
the corresponding row and column of `Omega`. Nearly singular LD matrices are
rejected; pruning or a scientifically justified regularization/PCA procedure
must occur before this API.

The MR layer does not model covariance caused by overlapping exposure and
outcome samples, weak-instrument bias, winner's curse, or uncertainty in an
estimated LD reference. Normal approximations are used for causal-estimate and
intercept tests. Results should be interpreted jointly across estimators and
diagnostics, not by selecting the smallest p-value.

## Meta-analysis and meta-regression

Study estimates are treated as independent normal observations with known,
strictly positive sampling variances. Fixed-effect weights are `1 / v[i]` and
random-effects weights are `1 / (v[i] + tau^2)`. The supported tau-squared
estimators are generalized DerSimonian-Laird, Paule-Mandel (solving
`Q_E(tau^2) = k - p`), and profile restricted maximum likelihood. Tau-squared
is constrained to be nonnegative.

The reported Cochran `Q` for an intercept-only model and residual `Q_E` for
meta-regression use fixed sampling-variance weights. I-squared is the
nonnegative conventional `100 * (Q - df) / Q`; H-squared is bounded below at
one. Moderator `Q_M` is the Wald chi-square test for all moderator coefficients,
excluding the intercept when present. Heterogeneity R-squared is the
nonnegative proportional reduction in tau-squared from the intercept-only fit.

Normal inference is the default. Student-t inference uses `k-p` denominator
DF. Knapp-Hartung multiplies the coefficient covariance by `Q_E/(k-p)` and
uses the same t DF; the modified form caps that multiplier below at one.
Prediction intervals are returned only for random-effects pooling. Correlated
effects, multilevel meta-analysis, robust/sandwich variance, publication-bias
diagnostics, and effect-size transformation are outside this first engine and
must not be inferred from these results.

## Reproducibility

Backend reductions can differ in order and rounding. Tests use the deterministic
JDistlib CPU backend. Production results record the requested policy, concrete
provider, device description, acceleration state, and whether automatic routing
was active. Cross-provider comparisons should use numerical tolerances rather
than bit equality.
