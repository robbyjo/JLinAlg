# Expected responses, scenarios, and marginal effects

`ModelPredictions` gives GLM and GEE fits one response-scale reporting API. It
uses the fitted coefficient covariance and analytic delta-method gradients. The
examples below use a probit GLM, but the same methods accept a `GeeResult`.

## Fit a binary probit model

The response is zero/one, and the intercept remains an explicit design column:

```java
double[] outcome = {0, 0, 0, 1, 0, 1, 1, 1};
double[][] design = {
    {1, -1.5, 0}, {1, -1.0, 1}, {1, -.5, 0}, {1, 0, 1},
    {1, .5, 0}, {1, 1, 1}, {1, 1.5, 0}, {1, 2, 1}
};

GlmResult fit = Glm.fit(outcome, design, GlmFamilies.probit());
```

`GlmFamilies.probit()` is an alias of `binomialProbit()`. Both use
`binomial(probit)`: fitted means are standard-normal CDF values, fixed
dispersion is one, and coefficient inference uses Wald z tests. Prior weights
retain the ordinary binomial contract: for grouped data, pass observed
proportions as the response and integer trial counts as prior weights.

The implementation evaluates the smaller normal tail and its logarithm from
the linear predictor. Thus likelihood, deviance, residuals, and Fisher weights
remain usable when a displayed probability has rounded to one. It does not
invent a positive information floor once the derivative itself is outside the
representable floating-point range.

## Expected responses

The same probit fit and response-scale estimands are available from a numeric
table through the source-build CLI:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar glm-predict `
  --input cohort.tsv `
  --response case `
  --predictors exposure,age,sex `
  --family probit `
  --estimand expected `
  --out predicted-probabilities.tsv
```

For standardized contrasts, reuse the observed rows and set one design column
to each scenario value:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar glm-predict `
  --input cohort.tsv --response case --predictors exposure,age,sex `
  --family probit --estimand difference `
  --scenario-column exposure --first 1 --second 0

java -jar build/cli/jlinalg-0.3.6.jar glm-predict `
  --input cohort.tsv --response case --predictors exposure,age,sex `
  --family probit --estimand ame --ame-column age
```

Use `--averaging at-average` only when the synthetic average-covariate
estimand is intended; population averaging is the default.

```java
double[][] people = {
    {1, -.75, 0},
    {1,  .25, 1},
    {1, 1.25, 0}
};

ExpectedResponse[] expected = ModelPredictions.expectedResponses(
    fit, GlmFamilies.probit(), people);

double probability = expected[0].mean().estimate();
double meanSe = expected[0].mean().standardError();
double lower = expected[0].mean().confidenceLower();
double upper = expected[0].mean().confidenceUpper();
```

Each row reports the linear predictor and its standard error as well as the
expected response, its delta-method standard error, and a response-scale
confidence interval. An overload accepts link-scale offsets and a confidence
level.

These are intervals for the fitted conditional mean. They are not prediction
intervals for future binary, count, Gaussian, random-effect, or otherwise noisy
outcomes. Residual/process variance is deliberately absent from this API so it
cannot be mistaken for mean uncertainty.

Prediction calls validate the fitted-model contract before calculating an
estimand. The supplied family must match the family recorded by the fit, and
the GLM or GEE fit must have converged. For a rank-deficient GLM, every requested
design row and an average-marginal-effect coordinate must be estimable under the
fitted row space; GEE predictions retain the rank contract established when the
GEE was fitted.

## Scenario differences and risk ratios

To standardize an exposure contrast over an observed population, construct two
design matrices that differ only in the exposure setting:

```java
double[][] exposed = new double[design.length][];
double[][] unexposed = new double[design.length][];
for (int row = 0; row < design.length; row++) {
    exposed[row] = design[row].clone();
    unexposed[row] = design[row].clone();
    exposed[row][2] = 1.0;
    unexposed[row][2] = 0.0;
}

ScenarioDifference difference = ModelPredictions.scenarioDifference(
    fit, GlmFamilies.probit(), exposed, unexposed);
RiskRatio ratio = ModelPredictions.riskRatio(
    fit, GlmFamilies.probit(), exposed, unexposed);
```

The defaults use `Averaging.POPULATION_AVERAGE`: transform every design row to
the response scale and then average. This is standardized or g-computation
semantics. It is not the same as `Averaging.AT_AVERAGE_COVARIATES`, which first
averages each design column and evaluates one synthetic covariate row:

```java
ScenarioDifference atAverage = ModelPredictions.scenarioDifference(
    fit, GlmFamilies.probit(), exposed, unexposed,
    null, null, Averaging.AT_AVERAGE_COVARIATES, 0.95);
```

The distinction matters for nonlinear links. In both cases the delta-method
variance forms one joint gradient from the two scenarios. The returned
`covarianceBetweenScenarios()` is therefore included in the difference rather
than incorrectly treating two predictions from the same fit as independent.

Risk-ratio inference uses the gradient of the log ratio, then transforms its
confidence limits back with `exp`. The reported `standardError()` is the
response-ratio delta-method standard error, while `logStandardError()` records
the inference scale. Both scenario means must be positive; use a scenario
difference when a Gaussian or identity-link mean is nonpositive.

## Average marginal effects

For a continuous covariate that enters the linear predictor directly, request
its zero-based design-column index:

```java
AverageMarginalEffect ame = ModelPredictions.averageMarginalEffect(
    fit, GlmFamilies.probit(), design, 1);
```

For coefficient `beta[j]`, the population-average estimand is the mean of
`inverseLink'(x_i beta) * beta[j]`. Its covariance gradient includes both the
direct derivative with respect to `beta[j]` and inverse-link curvature.
Built-in families provide analytic curvature; source-compatible custom
families use the documented numerical derivative default.

`AT_AVERAGE_COVARIATES` instead evaluates the derivative at the mean design
row. For indicator variables and multi-column encodings, prefer an explicit
scenario difference. For polynomial, spline, interaction, or other transformed
columns, this method differentiates with respect to the supplied design column;
the caller must apply the transformation's chain rule or use scenario changes.

## GEE and reference distributions

Every method above has the same overload for `GeeResult`. It uses the selected
GEE covariance matrix, including robust or finite-cluster variants. GLMs with
estimated dispersion use residual Student-t critical values; GEE fits configured
for cluster-t inference use their cluster denominator degrees of freedom;
fixed-dispersion GLMs and asymptotic GEE use normal critical values.

The frozen numerical audit uses base R `stats::glm(...,
family=binomial(link="probit"))` and explicit matrix delta-method calculations.
Regenerate its reference output with:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' `
  src/test/R/probit-prediction-reference.R
```

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [J. A. Nelder and R. W. M. Wedderburn (1972) — Generalized linear models](../CITATIONS.md#nelder-wedderburn-1972)
- [C. I. Bliss (1934) — The method of probits](../CITATIONS.md#bliss-1934) — [PMID: 17813446](https://pubmed.ncbi.nlm.nih.gov/17813446/)
- [Kung-Yee Liang and Scott L. Zeger (1986) — Longitudinal data analysis using generalized linear models](../CITATIONS.md#liang-zeger-1986)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
