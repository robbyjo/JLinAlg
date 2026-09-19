# Statistical audit fixes and validation

Validated locally on 2026-09-19 after review of commits `812025d`, `8910094`,
`6b1ede1`, `d8c7e9e`, `56fc927`, and `1242d93`. This change fixes all ten
reported findings; it does not extend the previously documented scientific scope.

## Follow-up audit from e89f5f5

The remaining 18 previously unaudited commits from `e89f5f5` through `10befe5c`
were reviewed at `1dfcabb`. Five further findings were corrected on 2026-09-19:

| Finding | Correction | Regression evidence |
|---|---|---|
| Phenotype units changed exported zero-information decisions | Normalize the information check by the projection scale: inverse residual variance for OLS and the retained projection infinity norm for REML | Export/read-back at response multipliers 1, 1e-12 and 1e12 preserves p-values and the expected score, effect and covariance scaling; constant and nuisance-confounded variants remain zero |
| Isolated zero-weight observations broke kernel inference | Skip exactly zero query influences before auxiliary local fitting | An isolated x=100 observation no longer breaks an interior query; changing its response leaves inference unchanged; insufficient support at nonzero influence still rejects |
| YAML emitted extra values for presence-only switches | Command-specific switch arity; true enables, false omits, CLI overrides retain precedence | Switch matrix, inherited overrides, valued TWAS/PWAS booleans and an actual configured meta-regression with and without an intercept |
| Repeated LD validation and products across molecular models | Prepared batch shared by marginal and joint CLI/API inference | Analytic correlated joint test, immutable prepared state, existing reference fixtures, invalid-LD and dependent-model rejection |
| Quadratic categorical-only kernel inference | Compute exact stratum mean and HC3 variance in linear work | Mixed ordered/unordered stratum matches SE = sqrt(sum squared residuals)/(m-1) |

The score degeneracy threshold is relative to the projection and genotype norm,
consistent with the REML scan's numerical-information check. This corrects unit
dependence; it does not promise inference for numerically unresolved directions.
The related-export scaling regression uses an identity covariance with scaled
optimizer bounds; the existing balanced random-intercept export test separately
checks a nontrivial related covariance against ANOVA REML.

Prepared molecular work costs O(n³ + kn²), plus O(k²n + k³) when joint inference
is requested. Dense LD validation remains cubic. Continuous/mixed kernel HC3
still uses auxiliary local fits; only categorical-only inference is linear.
Statistical assumptions and unsupported extensions are unchanged.

### Follow-up performance and reproduction

Java 25, CPU, same Windows host and seeded diagnostic workloads. Single probe
invocations before/after (including JVM warm-up effects) measured:

| Workload | Before | After |
|---|---:|---:|
| Joint molecular test, 1,000 variants and 30 models | 13.578787 s | 0.411732 s |
| Categorical kernel inference, 20,000 rows | 5.196117 s | 0.003189 s |

The molecular chi-square differed from the independent cached calculation by
7.11e-15; categorical mean differed from direct sums by 3.47e-18 and SE agreed
at double precision. These timings are local diagnostics, not JMH guarantees.
The checked-in benchmark uses one warm-up and the median of three runs:

| Workload | Median |
|---|---:|
| Joint molecular test, 500 variants and 30 models | 0.035664 s |
| Joint molecular test, 1,000 variants and 30 models | 0.386514 s |
| Categorical kernel inference, 20,000 rows | 0.000594 s |

```powershell
.\gradlew.bat --no-daemon check javadoc executableJar --console=plain
.\gradlew.bat --no-daemon benchmarkSummaryKernelAudit --console=plain
```

New regression classes are `GaussianScoreWriterTest`, `ProjectConfigurationFlagsTest`
and `PredictedOmicsTest`; `NonparametricExtensionsTest` adds both kernel regressions.
Existing `XwasWorkflowTest` and `RareScoreRelatedTest` check end-to-end integration.
The full validation completed with 907 tests: 904 passed, three optional CHOLMOD
checks skipped, no failures or errors. Javadoc, executable packaging, site links
and the 117-reference citation index also passed.

## Correctness gates

| Finding | Change | Independent or analytical regression gate |
|---|---|---|
| NB standard errors on exact fits | Recompute final unscaled weighted-design covariance by QR; the same residual-floor issue is also avoided for exact Gaussian fits before moderation | Base-R inverse information: SE 0.23990279696049335 and p 0.0038612085215580121 |
| Missing MI parameter uncertainty | Resample observed rows for each conditional update; PMM donors also come from that bootstrap | Binary pooled variance near 0.25/100 and continuous pooled variance near sample variance/100 |
| Incorrect multicategory probabilities | One baseline-category multinomial likelihood and joint information matrix | Imputed marginal proportions near 0.6/0.3/0.1 across 300 independent imputations |
| PMM row-order bias | Shuffle tied distances for each recipient before choosing the nearest donor pool | Forward/reversed observed values 1..100 both preserve mean 50.5 and use the full range |
| Infinite EWAS tail statistics | Request upper-tail normal quantiles directly | Two p=1e-20 probes give Z 10.416778225085585, p 2.078759535665973e-25; opposite signs give Z=0, p=1 |
| Inference from failed NB fits | NaN SE/statistic/df/p for nonconverged or clipped-boundary fits; CLI retains internal p=1 for the complete BH family but outputs NaN | All-zero feature and CLI family-size regression |
| ACAT tail precision/overflow | Cotangent small tails and scaling before reciprocal operations | Singleton identity through subnormal inputs, weighted result near 1e-300, 90-digit mpmath reference |
| Two-feature voom crash | Bound nearest-neighbor count by available features | Two-feature/four-sample fit completes with finite p-values |
| Voom repeated full sorting | Sort locations once; interpolate the local-constant tricube trend at at most 2,049 query quantiles | Direct brute-force precision weights agree within 0.2% on the seeded 1,500-feature workload |
| Quadratic SVA KDE | Linear binning, zero-padded FFT convolution, interpolation; spacing at most bandwidth/32 | Multimodal and isolated-tail KDE agrees within 0.05% relative error against direct sums |

The ACAT scalar reference is 0.001953404405770000045044061972 for equal-weight
p-values `[.02,.0004,.2,.1,.8]`, from 90-digit cotangents and the Cauchy survival
function. The previous 0.0019534044057701767 fixture inherited tangent
cancellation; its tolerance was not loosened to accommodate the fix.

## Bootstrap imputation contract

Bootstrap regression supplies conditional model uncertainty, followed by
response/donor sampling. Binary responses share the coherent multinomial
implementation with unordered categories; categorical predictors use dummy
columns. Weak ridge regularization keeps separated bootstrap fits finite;
nonconverged conditional fits stop the run. CLI metadata records this algorithm,
ridge, donor policy, iterations and seed. Seeded results are deterministic but
intentionally differ from the earlier point-estimate implementation.

This is a bootstrap MI method, not an exact Bayesian parameter sampler or a
numerical clone of R mice. Bootstrap logistic imputation is also an explicit
[mice method](https://amices.org/mice/reference/mice.impute.logreg.boot.html).
The tests establish the reported scalar uncertainty failures are corrected;
they do not certify every conditional specification. Valid Rubin inference
still needs appropriate MAR assumptions, compatible analysis/imputation models,
adequate observed information, and sufficiently mixed chains. Multilevel,
survey-weighted and MNAR imputation remain outside this interface.

## Performance

Java 25 on the Windows validation host, CPU, median of three diagnostic
invocations. These are reproducible macrobenchmarks, not isolated JMH results.

| Workload | Audit before | Same diagnostic probe after |
|---|---:|---:|
| Voom, 4,000 features x 12 samples | 8.876799 s | 0.054621 s |
| SVA local-FDR calculation, 16,000 features | 1.415963 s | 0.002974 s |

The checked-in `InferenceAuditBenchmark` adds one warm-up and times complete
workflows. It measured:

| Complete workflow | Features | Samples | Median |
|---|---:|---:|---:|
| Voom | 1,000 | 12 | 0.015632 s |
| Voom | 4,000 | 12 | 0.037073 s |
| Voom | 20,000 | 12 | 0.184462 s |
| SVA, five iterations and two factors | 2,000 | 32 | 0.039214 s |
| SVA, five iterations and two factors | 16,000 | 32 | 0.472283 s |
| AutoSVA, twenty iterations and two factors | 2,000 | 32 | 0.098587 s |
| AutoSVA, twenty iterations and two factors | 16,000 | 32 | 1.089894 s |

Voom interpolation uses exact trend values on sorted query quantiles and exact
evaluation for up to 128 features. SVA evaluates at most 256 features directly;
constant inputs take linear work. The FFT grid is capped at 131,072 intervals:
an excessively narrow bandwidth raises an explicit error rather than silently
coarsening the grid or changing the estimator. The approximate paths are checked
against direct calculations; neither claims universal identity with R lowess or
R's local-FDR smoothing spline. Other factor estimation and sample-Gram costs
remain; the local-FDR improvement does not make all SVA work linear in samples.

## Reproduction

```powershell
.\gradlew.bat check --console=plain
.\gradlew.bat benchmarkInferenceAudit --console=plain
```

Focused regression classes: `EmpiricalBayesDifferentialTest`, `VoomTrendTest`,
`MultipleImputationTest`, `RegionLevelEwasTest`, `AcatTestsTest`,
`GaussianKdeTest`, and `XwasInferenceCliTest`. Existing differential,
confounder, MR, regression and enrichment references remain in the full suite.

The audit's independent fixed-dispersion NB reference can be reproduced with R:

```r
X <- cbind(1, rep(0:1, each=3))
mu <- c(rep(20,3), rep(40,3))
alpha <- 0.04883002798420158
se <- sqrt(solve(crossprod(X, X * (mu/(1+alpha*mu))))[2,2])
c(se=se, p=2*pnorm(abs(log(2)/se), lower.tail=FALSE))
```
