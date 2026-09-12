# Rare-variant cohort summaries and meta-analysis

JLinAlg separates participant-level summary generation from the choice of
association test. `rare-score` exports one cohort's variant scores and their
covariances. `rare-meta` combines independent cohorts and performs the explicitly
selected single-variant, burden, SKAT, or SKAT-O tests. Cohort group-test p-values
alone are insufficient input for this score-based workflow.

The initial scope is quantitative traits, diploid additive biallelic variants,
and independent cohorts on a compatible phenotype scale. The built-in exporter
fits an unrelated-sample Gaussian model. RAREMETALWORKER/rvtests summaries from
appropriately adjusted related-sample quantitative-trait models can be imported;
this does not adjust for relatives or overlapping participants *between* cohorts.
Binary-trait rare-case calibration, conditional analysis, multiallelic covariance,
variable-threshold tests, and heterogeneous-effect kernel meta-analysis are future
extensions. These commands do not implement Raremetal2's `--useExact` method.

## Choose the question and test

| Test | Alternative and interpretation | Main output |
|---|---|---|
| `single` | Association of each variant | Score, variance, beta, SE, signed Z, p, log p, cohort directions |
| `burden` | Association of a specified weighted allele burden; often powerful when many effects have the same direction | Burden beta, SE, signed Z, p, log p, directions |
| `skat` | Set association allowing mixed effect directions and many non-associated variants | Quadratic Q, p, log p, numerical calibration |
| `skat-o` | Adaptive combination of burden and SKAT, accounting for the search over combinations | Minimum component p, adjusted p, rho/component diagnostics, calibration and simulation budget |

No test is uniformly most powerful. Choose the primary test, windows/masks, MAF
threshold, weights, and multiplicity correction before inspecting association
results. SKAT-O adjusts its internal rho search, not searches across genes,
overlapping windows, masks, or alternative tests. JLinAlg outputs unadjusted
across-feature p-values; apply the chosen multiple-testing procedure downstream.

SKAT and SKAT-O do not define a single signed window effect. They output no beta,
SE, or direction column. Select `--test burden,skat` to request both, and use the
separate burden file for the explicitly specified burden effect. Its confidence
interval is not a confidence interval for the SKAT association. Equal burden
can cancel opposing effects even when SKAT detects a strong association.

## A complete small example

Run from the repository root after `gradlew executableJar`. The executable uses
the current project version, presently `build/cli/jlinalg-0.3.5.jar`; these commands
are unreleased and are not in the previously published v0.3.5 asset. Example
cohorts have only six synthetic people, so `--maf 0.5` is used to retain variants
for illustration. This is not a realistic rare-variant power study.

```shell
java -jar build/cli/jlinalg-0.3.5.jar rare-score --vcf examples/rare-meta/cohort-a.vcf --pheno examples/rare-meta/cohort-a.tsv --id sample --response trait --genome-build GRCh38 --cov-window 100 --out build/rare-demo/a
java -jar build/cli/jlinalg-0.3.5.jar rare-score --vcf examples/rare-meta/cohort-b.vcf --pheno examples/rare-meta/cohort-b.tsv --id sample --response trait --genome-build GRCh38 --cov-window 100 --out build/rare-demo/b

java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --test single --out build/rare-demo/single
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test burden --weights equal --maf 0.5 --cohort-results --out build/rare-demo/burden
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test burden --weights mb --maf 0.5 --out build/rare-demo/weighted
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test skat --weights beta --maf 0.5 --out build/rare-demo/skat
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test skat-o --weights beta --maf 0.5 --simulations 100000 --seed 1234 --out build/rare-demo/skato
```

Each selected test creates `PREFIX.TEST.tsv`. `PREFIX.log` records settings,
cohort order, input paths, UTC start/end timestamps, status, and readable elapsed
time. Existing files are never overwritten. `--cohort-results` adds cohort rows
with `scope=A`, `scope=B`, etc., alongside `scope=meta`; otherwise only meta rows
are written. It uses the same pooled mask and weights for cohort comparisons.
Those contributing cohorts are not independent replication of their meta-result.

Single-variant output includes, for example, `feature_id=1:10:A:G`, the number of
informative cohorts, and `direction=++` when both scores favor the G allele.
`+`, `-`, `0`, and `?` mean positive, negative, zero, and unavailable information.
The burden file has `weights`, `direction`, `beta`, `se`, `z`, and `p_value`;
the SKAT file replaces signed-effect fields with `q` and `calibration`;
SKAT-O reports `minimum_component_p`, `adjusted_p`, and `component_rho_p`.

## Java score-summary example

The runnable [RareMetaExample.java](../examples/rare-meta/RareMetaExample.java)
performs single-variant, equal-weight burden, weighted burden, SKAT, and SKAT-O
meta-analysis using primitive arrays. With a JDK and the source-built JAR:

```shell
java --class-path build/cli/jlinalg-0.3.5.jar examples/rare-meta/RareMetaExample.java
```

Its two synthetic independent cohorts have scores for the same two variants.
Scores must already share allele orientation and trait units. Covariances here
are on the score-information scale, with one row-major matrix per cohort:

```java
import org.jlinalg.settest.*;

double[][] scores = {{2, -3}, {4, 1}};
double[][] covariance = {{4, 1, 1, 9}, {16, 2, 2, 4}};
var pooled = ScoreMetaAnalysis.pool(scores, covariance, 1);
if (pooled == null) throw new IllegalStateException("No eligible variants");
var state = pooled.state();
double[] u = state.scores(), v = state.information();
var first = SummarySetTests.singleVariant("v1", u[0], v[0]);
var burden = SummarySetTests.burden("GENE1", state, new double[]{1, 1});
var weightedBurden = SummarySetTests.burden("GENE1", state, new double[]{1, 2});
var skat = SummarySetTests.skat("GENE1", state, new double[]{1, 2});
var skato = SummarySetTests.skatO("GENE1", state, new double[]{1, 2},
    SetTestOptions.defaults());
```

The pooled scores are `(6, -2)` and covariance is `[[20,3],[3,13]]`. Equal
burden beta is `4/39` with SE `1/sqrt(39)`. The weights `(1,2)` are illustrative
linear coefficients, not automatically estimated MAF weights; select weights and
variant masks before testing. SKAT squares those coefficients in its kernel.
The full example prints both single-variant results and directions, burden
beta/SE/p, SKAT Q/p, and SKAT-O adjusted p with its simulation budget and seed.
The default SKAT-O calibration has the same resolution limits as the CLI.

Represent unavailable variants with `Double.NaN` scores; corresponding covariance
entries are ignored. The final `pool` argument is the minimum informative cohort
count per variant. `pooled.indices()` maps retained variants to input columns;
subset weights in that order if variants were filtered. The arrays contain
unweighted, nuisance-adjusted scores, not per-cohort burden/SKAT p-values.
Use `RareMetalStudy` for indexed file reads and covariance rescaling; `rare-meta`
also handles file-level allele alignment, masks, and MAF weight construction.

For gene/probe effects with beta and SE, see the
[omics Java and CLI meta-analysis example](vignettes/meta-analysis.md#omics-example-expression-effects-across-cohorts).

## Cohort preparation and interoperability

The manifest is TSV with `cohort`, `scores`, and optional `covariance` columns.
Paths are relative to the manifest. Covariance files are needed for multivariant
tests, not single-variant score meta-analysis. Each name must be unique and must
not be `meta`. Reusing the same score file as two cohorts is rejected.

`rare-score` accepts VCF/BGZF VCF/BCF, aligns sample IDs to the phenotype table,
uses complete phenotype/covariate cases, fits an intercept plus optional numeric
`--covariates age,sex,pc1`, and mean imputes missing dosages on that fixed sample.
Scores and information use original phenotype units and null residual variance
`RSS/(N-p)`. This estimator choice is recorded; it is not a claim that every
RAREMETALWORKER model or transformation produces identical summaries. No inverse
normal transformation is silently applied. HWE is not computed and is recorded
as `NA`; use `--hwe 0` for these outputs. For fractional DS inputs, use RAREMETAL's
`--dosage` option if consuming the exported files in that external program.

The exporter writes BGZF and tabix indices with the historical RMW column order.
Its program header identifies `JLinAlg (RareMetalWorker-compatible)` because
RAREMETAL 4.15.1 detects the column layout from that substring. Covariances are
divided by `AnalyzedSamples` on disk, as required by RMW/rvtests format, and are
multiplied back on import. Diagonal covariance is checked against `SQRT_V_STAT^2`.
The reader also accepts a `##CovarianceScale=score` header for JLinAlg-only,
unscaled covariance files; do not send that extension directly to RAREMETAL.

Per-variant scores and covariance alone do not encode all modeling decisions.
Agree on phenotype units/transforms, covariates, allele normalization, reference
build, and ancestry/relatedness adjustment before sharing summaries. Explicit
conflicting `GenomeBuild` headers are rejected; historical files without headers
rely on the user's `--genome-build` declaration. No liftover or strand inference
is performed. Exact REF/ALT swaps flip scores and both covariance axes.

Indexed files are recommended for regional access. Unindexed files are rescanned
per group with bounded memory, which is slower. Sequential scans require unique
positions in natural chromosome order (numeric autosomes, X, Y, MT, other contigs).
Ambiguous multiallelic positions and incompatible allele pairs fail. Historical
monomorphic rows with one known allele are interpreted explicitly; polymorphic
rows with an unknown alternate allele are excluded from single-variant inference.

## Windows, missingness, and weights

Use either a group file (`GROUP_ID CHROM:POS:REF:ALT ...`) or
`--window-size 10000 --window-step 5000`. Windows use 1-based inclusive positions,
start at 1 + a multiple of the step, and are evaluated if they contain variants.
Genes/groups must lie on one chromosome. The exporter must save covariance over
at least the maximum within-group distance. A required missing covariance pair
is an error, never an assumed zero. Explicit zero covariances are valid.

`--min-cohorts` defaults to 1 and counts informative cohorts **per variant**.
Groups retain variants meeting that threshold and the pooled MAF cutoff.
The result's `n_cohorts` counts cohorts contributing at least one retained variant;
it does not assert that every cohort measured every retained variant. No eligible
variants, unavailable information, and one-cohort pass-through have explicit
statuses. All-common, monomorphic, or QC-excluded groups are retained as status
rows rather than invented test results. Zero-information burdens get p=1 and
undefined beta/SE, with an explicit status.

`--af-policy observed` (default) pools available, QC-passing allele frequencies
weighted by informative sample counts; absent records contribute no denominator.
`--af-policy raremetal` additionally treats an absent record as reference-only in
that cohort for AF pooling, reproducing the historical convention on compatible
inputs. Absence is not proof of monomorphism; select this policy deliberately.
Group alleles are oriented to the pooled minor allele. Single-variant estimates
remain oriented to the first available complete REF/ALT pair.

Weights are linear coefficients `w` for burden `B=G w`; SKAT uses `w^2` on its
kernel diagonal. Options are `equal` (1), `mb` (1/sqrt(f(1-f))), `beta`
(Beta(1,25) density), and `raremetal-beta` (density with RAREMETAL's frequency
truncation at 2/total cohort N). Default burden weights are equal; default SKAT
and SKAT-O use `raremetal-beta`. These choices are explicit in output. `beta`
burden uses the density directly; it does not emulate RAREMETAL's unusual
inverse-squared-density `--BBeta` burden implementation.

For aligned independent cohorts, `U=sum(U_k)` and `V=sum(V_k)`. Single-variant
beta/SE are `U/V` and `1/sqrt(V)`; burden beta/SE are `(w'U)/(w'Vw)` and
`1/sqrt(w'Vw)`. These are one-step/null-information estimates, not universally
the alternative-model fitted coefficient and SE. Cohort burden estimates can
be used in ordinary fixed/random-effects meta-analysis only when the burden
definitions and units match; doing so is a separate model from kernel pooling.

## Calibration, numerical limits, and validation

Single/burden inference uses asymptotic normal scores and retains log p-values
even when ordinary double p-values underflow. SKAT uses deterministic positive
chi-square-mixture tails with the existing labeled saddlepoint fallback where
the bounded series cannot certify accuracy. Input covariance is checked for
symmetry, PSD, and score consistency with its numerical null space.

Default SKAT-O uses correlated Gaussian score-null simulation, in bounded blocks,
with an explicit seed and `(exceedances+1)/(simulations+1)` p-value. Its minimum
resolution is `1/(simulations+1)`; 10,000 simulations cannot resolve genome-wide
small p-values. Increase the budget for a prespecified precision target, or
explicitly request `--skato-calibration analytic`, a faster **moment approximation**.
Neither option establishes exact finite-sample phenotype calibration. Rank-one
SKAT-O reduces to a one-degree-of-freedom test without simulation.

See [validation and timing](rare-meta-validation.md) for the tested scope,
RAREMETAL version, independent R fixtures, and measured performance. Leave-one-out
and conditional analysis remain future diagnostics; neither is needed to recover
single-variant scores from these summary files.

## References

- [RAREMETAL 4.15.1 release](https://github.com/statgen/raremetal/releases/tag/v4.15.1)
- [RAREMETAL workflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC4173011/)
- [Score-based gene meta-analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC3939031/)
- [MetaSKAT framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC3710762/)
- [SKAT-O method](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
