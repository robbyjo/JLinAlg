# Regression inference validation and supported scope

Date: 2026-09-19. These source-build additions finish the five remaining
Zelig-inspired regression TODOs on the synchronized v0.3.6 code base.
The previously published probit, prediction/contrast, IV, xWAS and ARIMA
workflows retain their existing APIs. Eight research items remain open in
[TODO](../TODO.md); this change does not claim to complete that backlog.

## Numerical references

`src/test/resources/r-reference/generate-regression-inference.R` generates
the checked-in tables and `reference.properties` under
`src/test/resources/r-reference/regression-inference`. The deterministic seed
is 20260919; each fixture has 240 rows. Versions used: R 4.6.1, survival 3.8-6,
MASS 7.3-65, and survey 4.5. R is required only to regenerate references.
The generator prepends `build/r-library` for optional local packages.

| Method | Independent comparison | Absolute test tolerance |
| --- | --- | --- |
| Censored Gaussian | `survival::survreg`, left/exact/right observations | beta/scale 2e-6; covariance/log likelihood 2e-7 |
| Lognormal, Weibull, exponential AFT | `survival::survreg`, right censoring | beta/scale 2e-6; covariance/log likelihood 2e-7 |
| Ordered logit/probit | `MASS::polr` with fitted ordered thresholds | beta/threshold/covariance 3e-6; log likelihood 1e-7 |
| Rare-events logit | Independent R King–Zeng weighted-bias and prior-shift equations; ordinary ML comparison | coefficients/covariance 2e-7 |
| Survey-design regression | `survey::svyglm`, Gaussian/logit/probit/Poisson | coefficients/covariance 2e-6 |

These are fixture tolerances, not universal accuracy bounds. Zelig is a workflow
reference; the tests use underlying estimators or separately implemented equations.
The rare-event reference tests independently assembled correction equations,
not a claimed full Zelig execution. Additional checks cover singleton-stratum
rejection, survey weight-scale invariance, AFT time-unit equivariance,
ordered-probability normalization, separated-fit rejection, and invalid data.
CLI tests exercise every new command, prediction files, metadata, run logs,
invalid controls and output-collision preservation. All 11 executable commands
in the two tutorials were also run against the assembled source-build JAR.

The full `check executableJar` gate passed on 2026-09-19: 873 tests discovered,
870 passed, zero failures/errors, and three optional native CHOLMOD tests
skipped because that backend was unavailable. The eight new test methods all
passed. Website structure/link checks, citation-index consistency, Javadoc,
and the regression benchmark also passed.

Observed-information covariance for censored/ordinal ML differences the
**analytic score**, checks two step sizes, requires positive-definite information,
checks the inverse residual, and requires a small Newton step. It does not add a
ridge or substitute a pseudoinverse for unidentified uncertainty.

## Reproduce

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-regression-inference.R
New-Item -ItemType Directory -Force build/regression-temp | Out-Null
$env:JAVA_TOOL_OPTIONS='-Djava.io.tmpdir=E:/Projects/JLinAlg/build/regression-temp'
$env:GRADLE_USER_HOME='E:/Projects/JLinAlg/build/validation-gradle-home'
.\gradlew.bat --no-daemon check executableJar benchmarkRegressionInference
```

The workspace-specific temporary paths avoid Windows temporary-directory access
failures observed on this host. Other installations can use their ordinary
Gradle user home and JVM temporary directory. The usual portable gate is
`./gradlew check executableJar benchmarkRegressionInference`.

## Measured workload

CPU, Java 25 on Windows, 240 observations with two numeric predictors; fixtures
are already in memory. Three warmups and median of seven repetitions per method,
JVM heap 64 MiB initial / 1 GiB maximum. These measurements exclude CLI startup,
file parsing and reference generation. They do not establish large-data scaling
or a speed advantage over R.

| Operation | Median milliseconds |
| --- | ---: |
| Censored Gaussian | 0.7083 |
| Lognormal AFT | 0.6404 |
| Weibull AFT | 0.4347 |
| Exponential AFT | 0.3692 |
| Ordered logit | 0.9015 |
| Ordered probit | 1.3342 |
| Rare-events coefficient/prior correction | 2.3124 |
| Survey Gaussian | 0.6201 |

## Boundaries

- Censored models support left, exact and right observations, with independent
  censoring given covariates. No interval censoring, truncation, delayed entry,
  frailty, observation weights or robust covariance are supplied here.
- Ordinal fits assume common slopes and nonempty prespecified categories.
  Thresholds supply location, so the predictor design excludes an intercept.
- Rare-event correction requires a successful sample logistic ML fit. Fisher
  covariance is first-order, with supplied population prevalence treated as
  known. This is not finite-sample rare-variant tail calibration.
- Survey variance is with-replacement stratified PSU Taylor linearization, with
  no finite population correction. Singleton strata fail. Tests explicitly use
  PSUs minus strata df. Domain, replicate-weight and multistage inference remain
  outside this contract.

See [censored/ordinal tutorial](vignettes/censored-ordinal.md),
[sampling-model tutorial](vignettes/sampling-models.md), and
[method references](CITATIONS.md).
