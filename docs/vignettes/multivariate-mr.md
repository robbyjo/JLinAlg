# Multivariable and multivariate Mendelian randomization

The words *multivariable* and *multivariate* describe different models in
JLinAlg. The CLI therefore uses different commands and rejects a one-outcome
table at the multivariate boundary.

| Command | Exposures | Outcomes | Model |
|---|---:|---:|---|
| `mr-mvmr` | two or more | exactly one | multivariable MR (MVMR) direct effects |
| `mr-multivariate` | one or more | two or more, correlated | covariance-aware multivariate MR-IVW and multivariate MR-PRESSO |
| `mr-xwas` | many | many | bounded parallel **separate-pair** MR screening, not a joint model |

The new commands are available from a current source build. Build the
executable with `gradlew.bat assemble`; on this checkout the runnable artifact
is `build/cli/jlinalg-0.3.6.jar`. A published release JAR predating these
commands will not acquire them merely because the website has been updated.

## Prepare one wide, harmonized table

Clump the union of candidate instruments once, then align every reported
effect to the same effect allele before constructing this table. The joint
commands deliberately do not guess allele orientation.

```text
variant_id  beta_exposure_bmi  se_exposure_bmi  beta_exposure_ldl  se_exposure_ldl  beta_outcome_cad  se_outcome_cad  beta_outcome_stroke  se_outcome_stroke
rs101       0.08               0.02              0.03               0.01             0.025             0.010           0.018                0.012
```

Names passed to `--exposures` and `--outcomes` select the corresponding
`beta_exposure_NAME`, `se_exposure_NAME`, `beta_outcome_NAME`, and
`se_outcome_NAME` columns. CSV and TSV inputs are accepted. Every row must have
all selected associations; silently changing the instrument set by outcome
would change the joint model.

## Run MVMR: multiple exposures, one outcome

```powershell
java -jar build/cli/jlinalg-0.3.6.jar mr-mvmr `
  --input harmonized-wide.tsv `
  --exposures bmi,ldl `
  --outcome cad `
  --method ivw `
  --backend cpu `
  --output mvmr-cad.tsv
```

Use `--method egger` for a multivariable Egger fit. Each result row is a direct
effect conditional on the other listed exposures. `marginal_f` is the mean
squared marginal exposure z statistic; it must not be reported as a
conditional F statistic. The Java covariance-aware MVMR API can calculate
conditional strength when full within-variant cross-exposure sampling
covariances are available.

## Supply outcome correlation

`mr-multivariate` requires a headerless correlation matrix. Rows and columns
must follow `--outcomes` exactly:

```text
1.00  0.42
0.42  1.00
```

This is the correlation between the **estimated SNP-outcome associations**,
not automatically the phenotypic correlation. Estimate it from an appropriate
shared sample, cross-trait score covariance, or another justified source.
Record the source and sample overlap. The matrix must be finite, symmetric,
positive definite, and have a unit diagonal.

## Run covariance-aware multivariate MR-IVW

```powershell
java -jar build/cli/jlinalg-0.3.6.jar mr-multivariate `
  --input harmonized-wide.tsv `
  --exposures bmi,ldl `
  --outcomes cad,stroke `
  --outcome-correlation outcome-correlation.tsv `
  --method ivw-fixed `
  --backend cpu `
  --output multivariate-mr.tsv `
  --plot multivariate-mr-forest.svg
```

Use `ivw-random` for multiplicative-random standard errors. The estimator
stacks the outcome associations, builds each variant's covariance as
`diag(se) R diag(se)`, and fits the exposure-by-outcome coefficient matrix by
generalized least squares. Output includes:

- one coefficient test for every exposure-outcome cell;
- a joint test for each exposure across all outcomes;
- a joint test for each outcome across all exposures;
- an overall joint Wald test;
- multivariate Cochran Q, degrees of freedom, p-value, and dispersion; and
- the marginal instrument-strength summary for each exposure.

`--plot` writes a dependency-free SVG forest plot of every exposure-outcome
effect and confidence interval. For PRESSO it shows the corrected fit when one
is estimable, otherwise the raw fit.

Outcome correlation affects coefficient covariance and joint inference. A
loop of univariate fits cannot reproduce these joint tests.

## Run multivariate MR-PRESSO

```powershell
java -jar build/cli/jlinalg-0.3.6.jar mr-multivariate `
  --input harmonized-wide.tsv `
  --exposures bmi,ldl `
  --outcomes cad,stroke `
  --outcome-correlation outcome-correlation.tsv `
  --method presso `
  --simulations 10000 `
  --seed 20260916 `
  --outlier-threshold 0.05 `
  --output multivariate-presso-effects.tsv `
  --diagnostics multivariate-presso-instruments.tsv
```

The diagnostic uses leave-one-instrument-out residual vectors, Mahalanobis
distance under the supplied outcome covariance, a reproducible parametric
bootstrap, and Bonferroni-adjusted instrument tests. The effects file reports
raw and, when estimable, outlier-corrected fits. The diagnostic file records
the global statistic and empirical p-value, every instrument's distance and
p-values, the simulation count, and seed. Report the removed variants; do not
treat outlier removal as proof that the remaining instruments satisfy the
exclusion restriction.

## High-throughput xWAS workflow

Use `mr-xwas` to screen many exposure-outcome pairs efficiently. It keeps
bounded parallel work, deterministic ordering, complete-family BH FDR, and
hit-only diagnostics. It does not estimate a multivariate outcome covariance.
For a related phenotype family:

1. clump each exposure once;
2. screen all pairs with `mr-xwas`;
3. choose a prespecified related-outcome set and construct the complete wide
   table for retained instruments;
4. estimate or import the outcome-association correlation matrix; and
5. use `mr-multivariate` as the covariance-aware joint follow-up.

Do not select the joint phenotype set solely from the same unadjusted pairwise
p-values and then interpret its global p-value as confirmatory.

## Java API

```java
MultivariateMrResult fit =
    MultivariateMendelianRandomization.fit(
        instruments,
        List.of("bmi", "ldl"),
        List.of("cad", "stroke"),
        outcomeCorrelation,
        false,
        BackendPolicy.CPU);

MultivariateMrPressoResult presso =
    MultivariateMrPresso.analyze(
        instruments,
        List.of("bmi", "ldl"),
        List.of("cad", "stroke"),
        outcomeCorrelation,
        0.05, 10_000, 20260916L,
        BackendPolicy.CPU);
```

Coefficient storage is outcome-major, with exposures varying fastest. Prefer
`result.effect(exposureIndex, outcomeIndex)` when indexing manually.

## Method boundaries

JLinAlg's `mr-multivariate` implements covariance-aware summary-data
multivariate MR-IVW and the Mahalanobis/bootstrap multivariate MR-PRESSO
diagnostic. It does **not** label these as either of the following distinct
methods:

- MR², the sparse Bayesian Gaussian-copula model with MCMC exposure selection
  and learned residual graph; or
- MRMO, the individual-level two-stage mixed-response method for combinations
  such as continuous and binary outcomes.

Those methods have different estimands, input data, algorithms, and
convergence diagnostics. A binary GWAS log-odds association can be supplied to
the summary-data command on its documented scale, but that does not turn the
analysis into individual-level MRMO.

## Numerical validation

The fixed-effect implementation is checked against an independent base-R
stacked generalized-least-squares calculation. The frozen fixture compares all
coefficients, standard errors, the full coefficient covariance, multivariate
Q, and joint Wald statistics. Separate tests verify reduction to independent
MVMR fits when outcome correlation is the identity, covariance propagation,
random-effect inflation, deterministic PRESSO simulation, invalid-correlation
failures, and both CLI schemas.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Stephen Burgess and Simon G. Thompson (2015) — Multivariable Mendelian randomization: the use of pleiotropic genetic variants to estimate causal effects](../CITATIONS.md#burgess-mvmr-2015) — [PMID: 25632051](https://pubmed.ncbi.nlm.nih.gov/25632051/) · [PMCID: PMC4325677](https://pmc.ncbi.nlm.nih.gov/articles/PMC4325677/)
- [Eleanor Sanderson, Wesley Spiller, and Jack Bowden (2021) — Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization](../CITATIONS.md#sanderson-mvmr-diagnostics-2021) — [PMID: 34338327](https://pubmed.ncbi.nlm.nih.gov/34338327/) · [PMCID: PMC9479726](https://pmc.ncbi.nlm.nih.gov/articles/PMC9479726/)
- [Yuankai Zhang et al. (2026) — Multivariate Mendelian randomization for joint inferences of correlated outcomes](../CITATIONS.md#zhang-multivariate-mr-2026) — [PMID: 42207415](https://pubmed.ncbi.nlm.nih.gov/42207415/)
- [Verena Zuber et al. (2023) — Multi-response Mendelian randomization: Identification of shared and distinct exposures for multimorbidity and multiple related disease outcomes](../CITATIONS.md#zuber-mr2-2023) — [PMID: 37419091](https://pubmed.ncbi.nlm.nih.gov/37419091/) · [PMCID: PMC10357504](https://pmc.ncbi.nlm.nih.gov/articles/PMC10357504/)
- [Yangqing Deng et al. (2023) — Two-stage multivariate Mendelian randomization on multiple outcomes with mixed distributions](../CITATIONS.md#deng-mrmo-2023) — [PMID: 37338962](https://pubmed.ncbi.nlm.nih.gov/37338962/) · [PMCID: PMC10515454](https://pmc.ncbi.nlm.nih.gov/articles/PMC10515454/)
- [Marie Verbanck et al. (2018) — Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization](../CITATIONS.md#verbanck-mr-presso-2018) — [PMID: 29686387](https://pubmed.ncbi.nlm.nih.gov/29686387/) · [PMCID: PMC6083837](https://pmc.ncbi.nlm.nih.gov/articles/PMC6083837/)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
