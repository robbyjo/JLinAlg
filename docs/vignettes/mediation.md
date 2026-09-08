# Linear mediation analysis

`MediationAnalysis` implements a frequentist Gaussian mediation analysis from
three ordinary least-squares fits. It does not use bootstrap resampling,
Monte Carlo draws, or posterior sampling.

## Model

For treatment `X`, mediator `M`, outcome `Y`, and optional covariates `C`, the
component models are:

```text
M = iM + a X + alpha C + eM
Y = iY + c' X + b M + gamma C + eY
Y = iT + c X + delta C + eT
```

The reported effects are:

- `a`: treatment-to-mediator path;
- `b`: mediator-to-outcome path adjusted for treatment;
- `indirect = a * b`;
- `direct = c'`;
- `total = c`.

The indirect standard error is the first-order Sobel/product-of-coefficients
delta-method estimate

```text
sqrt(b^2 Var(a) + a^2 Var(b))
```

and its p-value and confidence interval use the asymptotic standard-normal
approximation. The direct, total, `a`, and `b` path intervals and tests retain
the usual OLS residual degrees of freedom.

## Example

Covariates should not include an intercept; the mediation fitter adds one to
each component model:

```java
import org.jlinalg.mediation.MediationAnalysis;
import org.jlinalg.mediation.MediationResult;

double[] treatment = {-2, -1, 0, 1, 2, 3};
double[] mediator = {0.2, 1.1, 2.0, 3.2, 3.9, 5.1};
double[] outcome = {1.1, 1.9, 3.3, 4.8, 6.0, 7.7};
double[][] covariates = {
    {0}, {1}, {0}, {1}, {0}, {1}
};

MediationResult result = MediationAnalysis.fit(
    outcome, treatment, mediator, covariates);

double indirect = result.indirectEffect().estimate();
double indirectP = result.indirectEffect().pValue();
double direct = result.directEffect().estimate();
double total = result.totalEffect().estimate();
```

Use `MediationEffect.standardError()`, `statistic()`, `pValue()`,
`confidenceLower()`, and `confidenceUpper()` for the reported uncertainty.
The component `OlsResult` objects remain available through
`mediatorModel()`, `outcomeModel()`, and `totalModel()` for diagnostics.

## Missing rows and scope

The mediation fitter applies one common complete-row decision to outcome,
treatment, mediator, and covariates. The default is
`MissingDataPolicy.ERROR`; set `MissingDataPolicy.OMIT` in `OlsOptions` to
remove incomplete rows. The retained original row indices are available from
`MediationResult.retainedRows()`.

This first implementation is deliberately limited to continuous Gaussian
responses and linear paths. Binary/count outcomes, nonlinear links, clustered
or robust covariance, and sampling-based intervals remain separate design
decisions rather than being silently approximated here.

## Mixed effects and pedigree terms

Use `MediationAnalysis.fitMixed` for ordinary independent random-effect terms.
The same terms are included in the mediator, outcome, and total-effect fits.
The sparse REML structure is prepared once and warm-started across the three
fits:

```java
RandomEffectTerm subject = RandomEffectTerm.randomIntercept(
    "subject", subjectIds);
MediationMixedResult mixed = MediationAnalysis.fitMixed(
    outcome, treatment, mediator, covariates,
    List.of(subject), RemlOptions.defaults(), 0.95,
    BackendPolicy.CPU);
```

Pedigree effects use the sparse additive relationship precision already used by
the pedigree mixed-model API. Ordinary random terms can be added alongside the
pedigree term:

```java
PedigreeRandomEffectTerm animal = PedigreeRandomEffectTerm.of(
    "animal", observedIndividualIds, pedigree);
MediationMixedResult pedigreeFit = MediationAnalysis.fitPedigree(
    outcome, treatment, mediator, covariates,
    List.of(animal), List.of(subject), RemlOptions.defaults(), 0.95,
    BackendPolicy.CPU);
```

Mixed and pedigree results expose the three underlying
`SparseLinearMixedModelResult` objects, including variance components,
convergence status, random effects, and fixed-effect covariance. These fits
require finite complete rows. Component path tests retain the selected mixed
model's denominator degrees of freedom and covariance; the indirect-effect
test remains an asymptotic normal Sobel test, not a Satterthwaite/Kenward–Roger
test of the product.

The Sobel calculation assumes zero cross-model covariance between the estimated
`a` and `b` paths. Separate mediator/outcome fits do not estimate correlated
cross-equation random effects or random-slope products. Consequently this API
is not a general multilevel causal mediation estimator. Causal interpretation
requires the relevant treatment/mediator confounding assumptions. In a mixed
fit, different fitted covariance weights also mean `total = direct + indirect`
need not hold, unlike ordinary same-sample OLS with these linear paths.

## Accuracy and speed comparison

Confidence intervals evaluate the small upper-tail probability directly,
including confidence levels immediately below one; forming `(1 + level) / 2`
can round to one and incorrectly produce infinite normal or mixed-path t
intervals. Independent R `qnorm`/`qt` references cover this boundary.

The deterministic comparator is `src/benchmark/r/mediation_benchmark.R` and
the matching Java task is:

```powershell
.\gradlew.bat `
  '-Djlinalg.benchmark.mediation.rows=2000' `
  '-Djlinalg.benchmark.mediation.groups=100' `
  benchmarkMediation
```

The OLS rows compare base R `lm` point estimates with JLinAlg's analytic
`a`, `b`, indirect, direct, and total effects. The optional `lme4` block
compares mixed-effect fixed effects. Run the R script with the same row and
group counts on a host with R and `lme4`; timings are workload- and
backend-specific, so a speedup claim must use the paired outputs rather than
the Java compile result alone.

The September 2026 audit freezes base-R and `lme4` references for all five path
effects and the Sobel SE. Maximum effect/SE differences are `1e-15` for OLS
and `2.3e-9` for the mixed fixture. A reciprocal-unit stress test also prevents
overflow caused by squaring a large path before multiplying a small variance:
the implementation now combines standard-error contributions with `hypot`.
Undefined zero/zero tests remain NaN; invalid variances and unrepresentable
effects are rejected explicitly.

See the [audit evidence](../../src/benchmark/resources/remaining-model-audit-benchmark/audit-results.md)
for the exact reference generator, reproducible isolated checks, and measured
warm timings. Those fixtures use first-order Sobel inference, not bootstrap
or posterior simulation from an R mediation package.
