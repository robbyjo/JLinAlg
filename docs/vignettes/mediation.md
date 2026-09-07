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
