# Censored, survival-time, and ordinal regression

These source-build commands use maximum likelihood, analytic scores, and
observed-information covariance. They require identifiable designs and a
resolved score/Newton step. Failed convergence, unresolved information, and
separation do not produce successful inference outputs. Build with
`./gradlew executableJar` or `.\gradlew.bat executableJar` on Windows.

## Censored Gaussian regression (Tobit)

```shell
java -jar build/cli/jlinalg-0.3.6.jar censored-regression --input src/test/resources/r-reference/regression-inference/censored.tsv --response tobit --censor cens --predictors c,z1 --distribution gaussian --predict src/test/resources/r-reference/regression-inference/censored.tsv --lower -0.4 --upper 2.5 --out build/inference/tobit.tsv
```

The synthetic table contains a normally distributed latent response censored
at -0.4 and 2.5. The response column stores the observed value or censoring
limit. Censor codes are **-1 left, 0 exact, 1 right**. A censored row contributes
a tail probability, not a normal density at the limit. Limits may differ by row
in fitting because each row supplies its actual limit.

The coefficient table includes an intercept, slopes, and `log(scale)`.
`scale` is the residual SD, not the variance. Prediction output separates:

- `latent_mean`: the uncensored Gaussian location x beta;
- `observed_mean`: E[clip(Y, lower, upper)];
- left and right censoring probabilities.

Prediction limits are fixed by `--lower`/`--upper`; an omitted side is unbounded.
Fitting does not infer censoring from equality with a limit: supply the codes.
Do not use this model for truncation or for a structural point mass at zero.

## Parametric accelerated-failure-time models

```shell
java -jar build/cli/jlinalg-0.3.6.jar censored-regression --input src/test/resources/r-reference/regression-inference/censored.tsv --response time --censor censTime --predictors c,z1 --distribution weibull --predict src/test/resources/r-reference/regression-inference/censored.tsv --time 2 --out build/inference/weibull.tsv
java -jar build/cli/jlinalg-0.3.6.jar censored-regression --input src/test/resources/r-reference/regression-inference/censored.tsv --response time --censor censTime --predictors c,z1 --distribution lognormal --out build/inference/lognormal.tsv
java -jar build/cli/jlinalg-0.3.6.jar censored-regression --input src/test/resources/r-reference/regression-inference/censored.tsv --response time --censor censTime --predictors c,z1 --distribution exponential --out build/inference/exponential.tsv
```

Times and censoring limits must be positive. Here `censTime=1` means right
censored, **not an event**; event indicators must be converted before use.
All models have `log(T) = x beta + scale * error`. Lognormal errors are standard
normal. Weibull errors have CDF `1-exp(-exp(error))`; therefore Weibull shape is
`1/scale` and its conventional scale is `exp(x beta)`. Exponential fixes the
log-time scale to one and estimates only regression coefficients.

Prediction output contains the median time, a 95% confidence interval for that
estimated median (including fitted scale uncertainty), and survival probability
at `--time`. The median confidence interval is not a future survival-time
prediction interval. A one-unit change in a direct predictor multiplies all
time quantiles by `exp(beta)`; it is a time ratio rather than a hazard ratio.

Java:

```java
CensoredRegression.Result fit = CensoredRegression.fit(time, design, censor,
    CensoredRegression.Distribution.WEIBULL);
double median = fit.timeQuantile(row, 0.5);
double survival = fit.survival(row, 2.0);
double[] medianCi = fit.timeQuantileInterval(row, 0.5, 0.95);
double[] timeRatio = fit.timeRatio(scenarioA, scenarioB, 0.95);
```

Covariance order is coefficients followed by log(scale), except exponential.
The Gaussian API additionally exposes `location`, `observedMean`, and
`leftCensoringProbability`. No interval censoring, delayed entry, truncation,
frailty, competing risks, robust sandwich covariance, or observation weights
are supplied by this new fitter. Censoring must be independent conditional on
the design. At least two exact observations are required; this does not alone
guarantee identification.

## Ordered logit or probit by maximum likelihood

```shell
java -jar build/cli/jlinalg-0.3.6.jar ordinal-regression --input src/test/resources/r-reference/regression-inference/ordinal.tsv --response ordinal --predictors c,z1 --categories 4 --link logit --predict src/test/resources/r-reference/regression-inference/ordinal.tsv --out build/inference/ordered-logit.tsv
java -jar build/cli/jlinalg-0.3.6.jar ordinal-regression --input src/test/resources/r-reference/regression-inference/ordinal.tsv --response ordinal --predictors c,z1 --categories 4 --link probit --out build/inference/ordered-probit.tsv
```

Responses are prespecified integers 0 through K-1, with every category observed.
The likelihood is `P(Y <= k | x) = F(threshold[k] - x beta)`. Increasing a
predictor with a positive coefficient shifts probability toward higher categories.
The link is logistic or normal. Threshold gaps use positive exponential
coordinates, and reported thresholds/covariance are transformed back to their
original coordinates. Numeric predictors **exclude an intercept**: thresholds
supply location. An intercept or constant in their span is rejected.

The model assumes common slopes across category boundaries (proportional odds
for logit). It differs from ordinal GEE, adjacent-category logits and ordinal
SEM. The `.predictions.tsv` reports each category's probability in row order.
Java uses `OrdinalRegression.fit(y, predictors, K, Link.LOGIT)` and
`fit.probabilities(row)`; covariance order is slopes then thresholds.

See [validation](../regression-inference-validation.md) for `survival::survreg`
and `MASS::polr` comparisons and [references](../CITATIONS.md) for methods.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [James Tobin (1958) — Estimation of Relationships for Limited Dependent Variables](../CITATIONS.md#tobin-1958)
- [John D. Kalbfleisch and Ross L. Prentice (2002) — The Statistical Analysis of Failure Time Data, Second Edition](../CITATIONS.md#kalbfleisch-prentice-2002)
- [Peter McCullagh (1980) — Regression Models for Ordinal Data](../CITATIONS.md#mccullagh-1980)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
