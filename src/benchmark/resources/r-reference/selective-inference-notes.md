# Fixed-penalty selective inference: validation and timing

The benchmark uses the existing deterministic n=160, p=4 regression fixture,
lambda=0.12, alpha=1, intercept included, no predictor standardization, known
noise SD=0.5, and 95% conditional confidence intervals. The active predictors
are columns 0 and 1. Both implementations target selected-model OLS
coefficients conditional on the selected active set and signs.

The selection objective is RSS/(2*n) + lambda*sum(abs(beta)). Accordingly,
`glmnet` and Java receive lambda=0.12, while `fixedLassoInf` receives
n*lambda=19.2. Lambda is fixed in advance and is not cross-validated.

## Reproduce

Run from the repository root:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/selective_inference_benchmark.R
javac -Xlint:all -Werror -cp 'build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar' -d build/selective-classes src/main/java/org/jlinalg/penalized/PolyhedralSelectiveInference.java
javac -Xlint:all -Werror -cp 'build/selective-classes;build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar' -d build/selective-benchmark-classes src/benchmark/java/org/jlinalg/benchmark/SelectiveInferenceBenchmark.java
java -cp 'build/selective-benchmark-classes;build/selective-classes;build/classes/java/main;build/dependencies/jdistlib-all-0.10.1.jar' org.jlinalg.benchmark.SelectiveInferenceBenchmark
```

The R script prefixes `.libPaths(c('build/r-library', .libPaths()))`. On the
validation host it needed elevated read access to the existing user-library
dependencies; no packages were installed. R 4.6.1, glmnet 5.0,
selectiveInference 1.2.5 were used. The script regenerates
`selective-n160.csv` beside this note.

`SELECTIVE_REPETITIONS` controls the R batch size; Java uses
`-Djlinalg.benchmark.selective.repetitions=...`. Both default to 200 fits per
batch and report the median of three batches. R uses three warmup fits per
mode; Java uses 100 warmup fits. Startup, package imports, compilation and
reference loading are excluded. Each timed fit includes penalized selection,
conditional inference and a checksum of beta, untruncated SE, both truncation
bounds, conditional p-value and both confidence endpoints for every effect.

## Accuracy gate and observed times

Before Java timing, active-set identity/count must match R. Each coefficient,
SE, truncation bound and tightly inverted CI endpoint must agree within 1e-8
absolute error. P-values, including 4.4051e-66 and 8.2211e-26, must agree
within 1e-8 relative error; an absolute tolerance that accepts zero is not used.
The R generator separately checks the selected-model OLS projection, KKT
normalization and attained CI tails. Both loops consume numeric output.

Windows host, 2026-09-08, 200 fits per batch:

| Workflow | Median seconds | Batch checksum |
| --- | ---: | ---: |
| Java, tight conditional-distribution inversion | 0.0249864 | -291.38904954668794 |
| R glmnet + native fixedLassoInf grid intervals | 0.400 | -291.16625461016349 |
| R native workflow followed by tight inversion | 0.510 | -291.38904954672472 |

The maximum absolute Java/reference difference was 2.26707542e-13 across
checked output fields; the separate relative p-value gates also passed.
These are one small deterministic workload's warmed timings, not universal
package speed ratios. The R tight workflow includes the native grid calculation
before replacing its intervals, so its time includes that additional work.

Native `fixedLassoInf` uses a CI search grid. Its endpoints differ from tight
inversion of the same truncation distribution by up to 0.00080546617901861017
on this fixture. The CSV preserves both sets (`grid_lower/grid_upper` versus
`ci_lower/ci_upper`). Java is checked against tight R `uniroot` inversion,
not against grid endpoints. "Tight" means numerical inversion with an
attained-tail check; neither implementation is claimed to use exact arithmetic.

## Boundary and scaling regressions

`SelectiveBoundaryTest` also checks cases absent from the n160 fixture:

- Five constant predictors, lambda=.4, and response just above the selection
  threshold (.4 + 1e-8, 1e-10, or 1e-12). Independent Mills-ratio/exponential
  tail limits verify the formerly erroneous confidence endpoints; sign
  reflection and attained conditional tail probabilities are checked.
- Finite and one-sided truncation with mean magnitude 1e12 and boundary gaps
  1e-12, plus intervals around zero narrower than ordinary CDF subtraction can
  resolve. Analytic truncated exponential and locally uniform limits apply.
- Extreme-tail sign symmetry with representable nonzero p-values; continuity
  and tiny-gap checks around the Mills-ratio branch transition.
- Dominating ridge penalties and predictor scales down to 1e-16. The formerly
  discarded lower bound is retained; p=0.13660884704067935 replaces the
  incorrect untruncated p=0.025347318677468245.
- Nonorthogonal n160 LASSO fits retain inference under predictor-unit changes
  by factors 1e-8 and 1e8, with the L1 penalty rescaled consistently.
- Unrepresentable variance fails explicitly, and confidence nearest one is
  inverted without rounding a target tail to exactly one.

The implementation normalizes inequalities, tests projection cancellation
relative to dot-product magnitude, and uses compensated dot products. Tail
ratios preserve boundary gaps before shifting by a large mean, using factored
normal exponents/Mills ratios and short-interval hazard integration. Confidence
endpoints require a valid bracket and an attained log-odds check; numerically
unresolved endpoints throw instead of being reported as successful intervals.

Known Gaussian noise variance, fixed penalty, the asserted selection event,
and identifiable selected-model projections remain assumptions. This work
does not provide unknown-variance or data-selected-lambda selective inference.
