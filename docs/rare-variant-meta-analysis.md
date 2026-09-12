# Rare-variant cohort summaries and meta-analysis

JLinAlg separates participant-level summary generation from the choice of
association test. `rare-score` exports one cohort's variant scores and their
covariances. `rare-meta` combines independent cohorts and performs the explicitly
selected single-variant, burden, SKAT, SKAT-O, variable-threshold, or heterogeneous-effect tests. Cohort group-test p-values
alone are insufficient input for this score-based workflow.

The scope is quantitative traits, diploid additive biallelic variants, and
independent cohorts on a compatible phenotype scale. The built-in exporter fits
an unrelated Gaussian model or a related-sample Gaussian REML model with `--grm`.
RAREMETALWORKER/rvtests quantitative-trait summaries can also be imported. This
does not adjust for relatives or overlapping participants **between** cohorts.
Binary-trait rare-case calibration, multiallelic covariance, and Raremetal2's
`--useExact` method remain outside the supported scope.
The [Raremetal2 assessment and trait-model contract](raremetal2-trait-models.md)
documents the required extensions and the implemented metadata compatibility checks.

## Choose the question and test

| Test | Alternative and interpretation | Main output |
|---|---|---|
| `single` | Association of each variant | Score, variance, beta, SE, signed Z, p, log p, cohort directions |
| `burden` | Association of a specified weighted allele burden; often powerful when many effects have the same direction | Burden beta, SE, signed Z, p, log p, directions |
| `skat` | Set association allowing mixed effect directions and many non-associated variants | Quadratic Q, p, log p, numerical calibration |
| `skat-o` | Adaptive combination of burden and SKAT, accounting for the search over combinations | Minimum component p, adjusted p, rho/component diagnostics, calibration and simulation budget |

| `vt` | Adaptive MAF-threshold burden search, calibrated under correlated Gaussian scores | Selected MAF and descriptive burden beta/SE/p; separate adjusted p, Monte Carlo SE and budget |
| `het-skat`, `het-skat-o` | Variant effects may differ between independent cohorts | Heterogeneous Q or adjusted omnibus p; separate number of cohort-by-variant effect dimensions |
| `burden-fixed`, `burden-random` | Common or normally distributed cohort burden effects with a comparable burden definition | Beta, SE, normal p; Cochran Q/p, I-squared, and REML tau-squared for the random model |

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

Each selected test creates `PREFIX.TEST.tsv`; `PREFIX.qc.tsv` records group coverage and exclusions. `PREFIX.log` records settings,
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
normal transformation is silently applied. HWE uses an exact two-sided test for
hard calls before imputation; fractional dosages or no calls produce `NA`.
A positive `--hwe` cutoff excludes rows with unavailable HWE. HWE in related or
ancestrally mixed samples is descriptive and its usual population interpretation
requires care; it is not adjusted by the GRM. For fractional DS inputs, use RAREMETAL's
`--dosage` option if consuming the exported files in that external program.

The exporter writes BGZF and tabix indices with the historical RMW column order.
Its program header identifies `JLinAlg (RareMetalWorker-compatible)` because
RAREMETAL 4.15.1 detects the column layout from that substring. Covariances are
divided by `AnalyzedSamples` on disk, as required by RMW/rvtests format, and are
multiplied back on import. Diagonal covariance is checked against `SQRT_V_STAT^2`.
The reader also accepts a `##CovarianceScale=score` header for JLinAlg-only,
unscaled covariance files; do not send that extension directly to RAREMETAL.

`rare-score` now writes version-1 quantitative model metadata into both files.
Use `--trait-id height --trait-units cm` to declare shared identity/units; defaults
are the response-column name and `original`. These options do not transform the
phenotype. `rare-meta --model-metadata strict` requires complete declarations in
each supplied file. Default `legacy` mode permits missing historical fields and
logs the assumptions. Both modes reject declared unsupported models/calibrations,
conflicting trait identities/units/transforms, mismatched score/covariance
metadata, and unsupported declared covariance column layouts. Matching metadata
cannot establish correct cohort provenance or independence.

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
RAREMETAL version, independent R fixtures, and measured performance.
Neither leave-out diagnostics nor conditioning is needed to recover the original
single-variant scores: they are already in the input files.

## Advanced tests, conditioning, and diagnostics

After generating the two cohort files in the complete example, these commands
run directly from the repository root. Use a fresh output prefix on each run.

```shell
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test burden,skat,skat-o,vt,het-skat,het-skat-o,burden-fixed,burden-random --weights equal --maf 0.5 --simulations 100000 --cohort-results --leave-variant-out --leave-cohort-out --threads 2 --cache-mb 16 --out build/rare-demo/advanced
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/conditional-groups.txt --condition examples/rare-meta/conditions.txt --test burden,skat-o --weights equal --maf 0.5 --skato-calibration deterministic --out build/rare-demo/conditional
java -jar build/cli/jlinalg-0.3.5.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test skat-o --weights equal --maf 0.5 --skato-calibration deterministic --out build/rare-demo/deterministic
java --class-path build/cli/jlinalg-0.3.5.jar examples/rare-meta/AdvancedRareMetaExample.java
```

The [Java example](../examples/rare-meta/AdvancedRareMetaExample.java) demonstrates
all model building blocks, including `SummaryScoreModels.condition`, `subset`,
`heterogeneous`, and `variableThreshold`, plus existing `MetaAnalysis` APIs.

### Conditional summary analysis

Both group and condition files have `GROUP_ID CHROM:POS:REF:ALT ...` lines.
The example tests variant 1 given variant 2; a singleton target group therefore
provides a conditional single-variant score test through `--test burden` with
equal weights. Target and conditioning positions must be disjoint and lie on
one chromosome. Groups absent from the condition file remain unconditional;
unknown condition group IDs are errors. The single-site streaming command
rejects `--condition` to avoid implying that its unconditional output is adjusted.

Within **each cohort**, compute `Uc = Ut - Vtc Vcc^-1 Uc0` and
`Vc = Vtt - Vtc Vcc^-1 Vct`, then pool the adjusted blocks. The implementation
uses a symmetric eigensolve, rejects conditioning eigenvalues at or below
`1e-10` of the largest, and validates the joint and adjusted covariance.
These are one-step conditional score estimates based on the supplied null model,
not an alternative-model participant-level refit. Omission uses a principal
submatrix and is a different operation.

All target/condition and condition/condition covariance pairs must be available;
extend the export covariance window accordingly. The default missing-conditioning
policy is an error. `--condition-missing exclude` excludes that entire cohort
from the group's scores and reports it in the QC manifest. The initial MAF/weight
definition still uses the original available cohort AFs; the informative cohort
minimum excludes unavailable conditioning cohorts. Condition variants are not
subject to the target rare-MAF filter, but must pass cohort QC and have information.

### Leave-out influence diagnostics

`--leave-variant-out` and `--leave-cohort-out` add rows to each requested group-test
table with `scope=leave_variant:CHROM:POS:REF:ALT` or `scope=leave_cohort:NAME`.
The original MAFs, weights, total-N weight truncation, and complete-case cohort
null models remain fixed. The original minimum-cohort filter defines the mask;
it is not reapplied after omission. Variants losing all information are dropped
and the result reports the remaining variant/cohort counts. Removing the only
variant or informative cohort produces an explicit unavailable-result row.
Every reduced SKAT-O or VT search is recalibrated using its reduced covariance.
These are exploratory influence checks, not conditional or causal attribution.
They are available for group tests, including singleton groups, rather than the
streaming `single` command.

### Variable thresholds and heterogeneous effects

VT searches the distinct pooled MAFs within the retained mask up to `--maf`, using
the chosen fixed linear weights (equal by default). It maximizes absolute burden Z
and uses joint Gaussian score draws to calibrate the correlated threshold search.
The selected burden beta/SE/p are descriptive after selection; use `adjusted_p`
for the threshold search and correct across features separately. `mc_se` reports
simulation uncertainty and the minimum resolution is `1/(simulations+1)`.
One informative threshold reduces to an ordinary normal burden test.

Homogeneous SKAT uses `sum_j (w_j sum_c U_cj)^2`. Heterogeneous SKAT instead uses
`sum_c sum_j (w_j U_cj)^2`, with block-diagonal covariance across independent
cohorts. Heterogeneous SKAT-O combines that kernel with the **pooled burden** and
calibrates the rho search. Thus opposite cohort effects can cancel in the pooled
burden while remaining visible to the heterogeneous kernel. `n_variants` counts
distinct retained variants; `n_effect_dimensions` counts informative cohort-by-
variant entries. These kernel models do not estimate a signed window beta or
shrinkage effects from their p-values.

`burden-fixed` and `burden-random` use cohort burden beta/SE only for cohorts with
the complete original mask and positive burden information. Incomplete-mask
cohorts are excluded from these scalar fits, marked `?` in their direction string,
and identifiable through the QC manifest. Random effects use REML tau-squared
and normal inference; with one eligible cohort, heterogeneity is unavailable and
the cohort estimate passes through with `single_cohort` status. These choices
do not make different phenotype units comparable. The Java `MetaAnalysis` and
`MetaRegression` APIs support additional explicitly chosen inference/model options.

### Deterministic rare-tail SKAT-O

`--skato-calibration deterministic` conditions on the Gaussian burden coordinate
and integrates the resulting **noncentral** quadratic-form tails. It avoids the
moment approximation and finite simulation resolution. Positive gamma-series
truncation is bounded by `1e-10 * minimum_component_p`; scaled adaptive quadrature
targets absolute error `1e-6 * minimum_component_p`. Quadrature error is an
estimate, not an interval-arithmetic guarantee. These limits are fixed and
identified in the calibration label. Supported min-p is at least `1e-250`.

An unresolved spectrum, series, quadrature, or component tail fails explicitly;
this mode never silently substitutes the saddlepoint or moment approximation.
Ill-conditioned or large noncentral problems may exceed its 8,192-term budget.
Independent rank-two polar and rank-three convolution references validate tested
moderate and rare tails; this does not establish universal numerical accuracy or
exact finite-sample calibration for sparse or non-Gaussian phenotypes.

### Related-sample export, QC, and memory

```shell
java -jar build/cli/jlinalg-0.3.5.jar rare-score --vcf examples/rare-meta/cohort-a.vcf --pheno examples/rare-meta/cohort-a.tsv --id sample --response trait --genome-build GRCh38 --grm examples/rare-meta/cohort-a.grm.tsv --cov-window 100 --out build/rare-demo/a-related
```

`--grm` accepts a labeled dense matrix or GCTA binary prefix, as described in
the [GRM tutorial](grm-cli.md). IDs are aligned to the fixed complete-case
VCF/phenotype sample. The model is `Var(y)=tau*K + sigma2*I`, fitted by REML;
nonconvergence is an error. Exported `U=G'Py`, `V=G'PG` retain original trait
units and disk covariance remains `V/N`. The log records fitted variance components.
Missing dosages are mean imputed and HWE is computed before imputation. Dense
REML retains an N-by-N projection, so this exporter is not a sparse biobank solver.

`PREFIX.qc.tsv` records requested target and conditioning variants per cohort,
retention, exclusion reason, sample count, pooled MAF, original information,
effect orientation, conditioning status, call rate, and HWE p. It does not assert
that overlapping groups or their samples are independent. For conditioned tests,
the manifest's information is the **original** cohort diagonal; test results use
the adjusted block. A single-site-only run produces a header-only group manifest.

`--threads` defaults to 1 and caps concurrent groups; output stays in input order.
Shared file readers serialize regional I/O while numerical tests run concurrently.
`--cache-mb` defaults to 16 per cohort, bounds estimated cached payload bytes, and
0 disables it. This is not a total heap limit. Active dense blocks require
approximately O(threads * cohorts * variants^2) storage plus eigensolver workspaces;
the heterogeneous kernel dimension is the sum of cohort-specific variant counts.
`--max-variants` caps target-plus-condition and heterogeneous dimensions. The
group output buffer is bounded by concurrent groups, but requesting all leave-out
rows adds memory and substantial computation. Indexed BGZF/tabix inputs are
recommended; unindexed files are rescanned on cache misses.

## References

- [RAREMETAL 4.15.1 release](https://github.com/statgen/raremetal/releases/tag/v4.15.1)
- [RAREMETAL workflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC4173011/)
- [Score-based gene meta-analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC3939031/)
- [MetaSKAT framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC3710762/)
- [SKAT-O method](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
