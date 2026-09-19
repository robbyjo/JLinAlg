# Rare-event logistic and survey-design inference

Build these source additions with `./gradlew executableJar` (Windows:
`.\gradlew.bat executableJar`). The examples use a synthetic 240-row fixture.
Input is CSV/TSV with explicit numeric predictors and complete finite outcomes.
Both commands include an intercept and write coefficient, metadata and log files.

## Rare-event logistic coefficient and prevalence corrections

```shell
java -jar build/cli/jlinalg-0.3.6.jar rare-events-logit --input src/test/resources/r-reference/regression-inference/survey.tsv --response rare --predictors c,z1 --out build/inference/rare-bias.tsv
java -jar build/cli/jlinalg-0.3.6.jar rare-events-logit --input src/test/resources/r-reference/regression-inference/survey.tsv --response rare --predictors c,z1 --bias-correction false --prevalence 0.04 --out build/inference/rare-prior.tsv
java -jar build/cli/jlinalg-0.3.6.jar rare-events-logit --input src/test/resources/r-reference/regression-inference/survey.tsv --response rare --predictors c,z1 --prevalence 0.04 --out build/inference/rare-both.tsv
```

Responses must be individual 0/1 observations. `--bias-correction true` is the
default King–Zeng first-order coefficient correction. `--prevalence` separately
specifies a known population event proportion for prior correction under
case-control sampling. The intercept receives logit(population prevalence)
minus logit(sample prevalence). If both corrections are requested, the bias
equations incorporate case/control sampling weights before the prior shift.
They do not fit a weighted pseudo-likelihood. Metadata records the sample
prevalence, bias option and prior intercept shift.

SEs are the underlying sample ML Fisher covariance, a first-order approximation.
They do not propagate uncertainty in supplied prevalence or give a finite-sample
coverage guarantee after bias correction. This is not iterative Firth fitting,
a remedy for separation, or rare-variant group-tail calibration. The underlying
logistic fit must converge at full rank. Both outcome classes are required.

Java uses `RareEventsLogistic.fit(y, design, biasCorrection, prevalence)`.
Supply a unit intercept as the first design column and `null` to omit prevalence
adjustment. The result exposes the coefficient-bias vector and prior shift
separately from the corrected coefficients.

## Survey regression with weights, strata and PSUs

```shell
java -jar build/cli/jlinalg-0.3.6.jar survey-regression --input src/test/resources/r-reference/regression-inference/survey.tsv --response y --predictors c,z1 --weights weight --strata stratum --psu cluster --family gaussian --out build/inference/survey-gaussian.tsv
java -jar build/cli/jlinalg-0.3.6.jar survey-regression --input src/test/resources/r-reference/regression-inference/survey.tsv --response rare --predictors c,z1 --weights weight --strata stratum --psu cluster --family logit --out build/inference/survey-logit.tsv
```

Here `weight` contains positive sampling weights; `stratum` identifies sampling
strata and `cluster` identifies primary sampling units (PSUs). PSU labels are
nested within strata. Each stratum must contain at least two observed PSUs;
single-PSU strata fail rather than silently contributing zero variance.

The coefficient fit uses sampling weights rescaled to mean one. Covariance
aggregates linearized score contributions within PSU, centers those totals
within each stratum, and applies the stratum's m/(m-1) multiplier. Multiplying
every weight by a positive constant does not change estimates or covariance.

Supported families: `gaussian` (default), `logit`, `probit`, and `poisson`.
Binomial responses are individual 0/1 outcomes; sampling weights are not trial
counts. Coefficient tests and intervals use **PSUs minus strata** t degrees of
freedom, explicitly rather than `svyglm`'s default model-adjusted residual df.
The metadata records the df and design counts.

This is a with-replacement, one-stage Taylor design variance or an ultimate-PSU
approximation. It does not implement finite population corrections, certainty
PSUs, multistage variance, replicate weights, calibration/raking, or domain
estimation by subsetting. Ordinary `--weights` in another regression command
does not provide this design variance.

Java:

```java
var sampling = new SurveyRegression.Design(weights, strata, psu);
var fit = SurveyRegression.fit(y, design, GlmFamilies.gaussian(), sampling);
double[] covariance = fit.covariance();
int df = fit.degreesOfFreedom();
```

See [validation](../regression-inference-validation.md) and
[method references](../CITATIONS.md).

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Gary King and Langche Zeng (2001) — Logistic Regression in Rare Events Data](../CITATIONS.md#king-zeng-2001)
- [Gary King and Langche Zeng (2001) — Explaining Rare Events in International Relations](../CITATIONS.md#king-zeng-isq-2001)
- [Thomas Lumley (2004) — Analysis of Complex Survey Samples](../CITATIONS.md#lumley-survey-2004)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
