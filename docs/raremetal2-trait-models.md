# Raremetal2 assessment and trait-model contract

Assessed on 2026-09-12 against upstream commit
[`da432cb8a3a1170022eac9d4ba29514ee601abe2`](https://github.com/statgen/Raremetal2/tree/da432cb8a3a1170022eac9d4ba29514ee601abe2).
This completes the format/method assessment and metadata prerequisite in the
original TODO. It does **not** implement Raremetal2 parity or add a binary,
Poisson, or survival model to `rare-meta`.

The delivered code validates quantitative model declarations, checks cohort
compatibility, exports a versioned contract, and rejects unsupported covariance
headers. Numerical quantitative score generation and pooling are unchanged.

## Three separate extensions

| Extension | Information needed | Decision |
|---|---|---|
| Multiallelic groups | Allele identity on both covariance axes; full same-site and cross-site covariance; producer-specific compression and ordering | Keep the biallelic boundary until a dedicated, independently validated decoder and allele-aware union/group implementation exist. |
| Raremetal2 `--useExact` | Cohort phenotype means/variances, counts, allele frequencies, compatible residualization and covariance scaling | Keep separate from ordinary sum-of-score pooling and from rare-case calibration. Define and validate the intended pooled model first. |
| Binary/other traits | Model-specific null fit, effect scale, sampling design, and enough information to calibrate tails | Version 1 remains Gaussian. Design future model versions around calibrated distributions, not just a renamed `U,V` file. |

### Multiallelic and compressed covariance

The upstream [group reader](https://github.com/statgen/Raremetal2/blob/da432cb8a3a1170022eac9d4ba29514ee601abe2/raremetal/src/annoGroups.cpp)
recognizes the old `CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES` layout and a
six-column layout with `REF ALT EXP COV_MATRICES` after chromosome and position.
The latter associates alleles with each row, decodes comma-separated values or
`value:index` entries, uses row-order differences to identify pairs, and rescales
with powers of ten and sample size. A bounded regional decoder consequently
needs the global row ordering and matrix coverage, including rows excluded by
group membership or QC. Allele swapping must change both covariance axes.

The same checkout's [worker](https://github.com/statgen/Raremetal2/blob/da432cb8a3a1170022eac9d4ba29514ee601abe2/raremetalworker/src/PreMeta.cpp)
also contains a four-column compressed `CHROM CURRENT_POS EXP COV_MATRICES`
writer. This differs from the allele-bearing header expected by the group
reader. The presence of four or six fields alone is therefore insufficient to
identify RMW versus rvtests versus Raremetal2 encoding. These are source-level
observations, not a successful native producer/consumer interoperability test.

JLinAlg now checks declared covariance column names before reading a matrix.
It accepts the old RMW layout and rvtests
`CHROM START_POS END_POS NUM_MARKER MARKER_POS COV`. Declared `EXP`, allele-aware,
and unknown layouts fail with a dedicated-decoder diagnostic, including when
only single-site results were requested but a covariance file was supplied.
Headerless historical files retain the old four/six-column interpretation in
legacy mode; strict mode requires a covariance column header.

Repeated score positions still fail in sequential and regional access. An old
position-only covariance file cannot distinguish `1:20:C:T` from `1:20:C:G` or
their covariance. Neither selecting the first row nor filling unknown pairs
with zero is valid. Before enabling this extension, fixtures must cover two
ALTs at a site, cross-ALT covariance, swapped alleles, missing alleles, explicit
zeros versus uncovered pairs, dense/sparse compression, reordered group members,
and equality of indexed/unindexed and serial/cached/parallel results.

### What `--useExact` means here

The pinned build selects `Meta.cpp` through its
[Makefile](https://github.com/statgen/Raremetal2/blob/da432cb8a3a1170022eac9d4ba29514ee601abe2/raremetal/src/Makefile).
In [that implementation](https://github.com/statgen/Raremetal2/blob/da432cb8a3a1170022eac9d4ba29514ee601abe2/raremetal/src/Meta.cpp),
`updateYstat` reads `TraitSummaries`; setup calculates cohort deviations from a
sample-size-weighted phenotype mean. `adjustStatsForExact` changes pooled scores
and variances using those deviations, per-cohort frequencies, counts and
phenotype variances. `updateExactCov` is a separate group-covariance path.
Optional population normalization additionally uses reference frequencies.

This is not obtained by changing a Gaussian tail routine, combining cohort
p-values, or applying the current `--af-policy`. Our files contain residual
variance (and REML metadata), but do not contain the complete, defined phenotype
moment/residualization contract that this method needs. Raw outcome variance,
residual variance with `N-p` degrees of freedom, and related-sample variance
components are different quantities. No conversion among them is inferred.

An implementation needs a written estimand for the pooled model and separate
participant-level reference fits with unequal sample sizes, phenotype means and
variances, shifted allele frequencies, missing variants, covariates and allele
reversals. Group covariance must match independent cross-products, with rank
and PSD checks. Related samples and population normalization need separate
validation. The upstream repository describes this as beta software; a native
comparison would be useful evidence but would not replace those checks.

## Implemented version-1 metadata

Both `rare-score` files now carry the following declarations. They describe
the fitted Gaussian score model; they do not certify the scientific provenance
or independence of the cohorts.

| Header | Supported value or meaning |
|---|---|
| `SummaryMetadataVersion` | `1` |
| `TraitType` | `quantitative` |
| `TraitId` | Shared phenotype identity, defaulting to `--response`; override with `--trait-id`. |
| `TraitUnits` | Shared units, default `original`; override with `--trait-units`, for example `cm`. No conversion is performed. |
| `TraitTransformation` | Exported as `none`; imported nonblank labels must agree across cohorts. Historical `InverseNormal=ON/OFF` maps to `inverse-normal/none`. |
| `NullModel` | `unrelated-Gaussian` or `related-Gaussian-REML` |
| `GenotypeModel` | `diploid-additive` |
| `EffectScale` | `trait-per-alt-allele` |
| `ScoreCalibration` | `asymptotic-normal` |
| `CovarianceModel` | `model-based` |
| `CovarianceScale` | Exported as `per-sample` (`V/N`); importer also accepts `score` (`V`). |

`AnalyzedSamples` remains mandatory in the score file. `GenomeBuild` is checked
against the CLI declaration when present. Supported keys cannot be blank or
duplicated. Unknown producer headers are ignored. Declared unsupported trait,
null-model, genotype, effect, covariance, calibration or schema values always
fail. A `BinaryTrait=True` declaration also fails; historical producers do not
consistently emit such a flag, so its absence is not proof of a Gaussian fit.

Score/covariance declarations must agree within a cohort, including sample
count, build, null model and storage scale. A missing score-file declaration can
be supplied by its covariance file in legacy mode and participates in cohort
checks. Across cohorts, declared identities, units, transformations and score
semantics must agree. OLS and Gaussian REML may coexist on the same effect scale;
their residual variances need not match. Storage scale can differ between
cohorts because the reader restores each matrix to score scale first.

Use `rare-meta --model-metadata strict` to require all version-1 fields in each
supplied score and covariance file. The default `legacy` mode preserves old
quantitative imports and logs `legacy-assumed-quantitative` for incomplete
effective declarations. Complete declarations log `declared-quantitative`.
The log records the policy and each cohort's recognized declarations. Legacy
mode does not override conflicts or make unsupported models acceptable. Missing
fields cannot establish compatibility: users must inspect producer commands and
phenotype definitions. The label `original` is a user assertion of matching
units, not an automatic check that different laboratories used the same units.

For example, add `--trait-id height --trait-units cm` to each existing
`rare-score` command, then add `--model-metadata strict` to `rare-meta`.
Matching different response-column names requires an explicit shared trait ID.
The older Java `GaussianScoreWriter.write` overloads remain source compatible
but omit `TraitId`; use the new trait-ID/unit overloads for strict imports.

This contract does not yet encode a cryptographic pairing of score/covariance
files, a shared covariate design, actual outcome distributions, or participant
overlap. Matching headers alone cannot establish any of those properties.

## Binary rare-case calibration prerequisite

For an unrelated Bernoulli logistic null with fitted probabilities `mu`, let
`W=diag(mu*(1-mu))`. With null design `X` and dosage matrix `G`, efficient scores
use `U=G'(y-mu)` and information
`V=G'WG-G'WX(X'WX)^(-1)X'WG`, subject to null convergence and rank checks.
These two moments support an asymptotic normal approximation. They do not
describe the discrete rare-case tail.

A concrete example has 1,000 participants, one case, and ten heterozygous
carriers, including that case. The intercept-only logistic score is `U=.99`
and information is `.0098901`. Its normal tail is below `1e-20`, whereas
conditioning on the total case count gives probability `.01` for an absolute
score at least as extreme: the single case lands among ten of 1,000 carriers.
The two distributions differ by construction. This is an illustration of the
calibration boundary, not a proposal to apply a hypergeometric tail to adjusted
or related-sample logistic scores.

The deterministic [base-R generator](../src/benchmark/r/rare_binary_calibration_boundary.R)
and [recorded fixture](../src/test/resources/raremetal/binary-calibration-boundary.tsv)
compute both quantities independently, without sampling. R's
[hypergeometric distribution documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html)
specifies the conditional mass calculation.

A future binary schema must specify case/control mapping and counts, sampling
or ascertainment design, offsets and trial weights, covariates/conditioning set,
null-fit convergence and rank, effect units (`log-odds-per-alt-allele`), dosage
handling, relatedness, covariance coverage, and calibration method/version/status.
Failures, boundary outcomes, and separation must suppress calibrated inference.
Per-variant calibrated p-values do not calibrate a multivariant SKAT statistic.

For a fixed fitted independent Bernoulli null, a burden's cumulant generator
can be built from residualized dosage values and fitted probabilities.
Nuisance estimation, conditioning on case count, mixed effects and case-control
ascertainment change the calibration problem. Aggregate counts and `U,V` alone
do not retain that information. A cohort-side evaluator or a carefully scoped
richer summary protocol is needed; do not export participant probabilities or
dosages under an aggregate-only interface. Independent-cohort cumulants add
only under a specified compatible model. Joint group tails and overlap require
additional information and validation.

Acceptance gates should include exact finite-support cases, independent
logistic likelihood/score/information checks, sparse-carrier and unbalanced
case-control fixtures, extreme-tail stability, explicit unattainable-tail and
nonconvergence behavior, and controlled calibration experiments across prespecified
case/MAC/covariate regimes. Numerical integration accuracy and empirical
type-I-error calibration must be reported separately.

Poisson models additionally require exposure/offset and dispersion semantics;
Cox models require censoring/truncation, ties, strata and risk-set information.
The existing [conditional-GWAS exports](conditional-gwas-summary.md) already
provide model-specific normal-only scores in a **different schema**. They are
not imports for `rare-meta`, and their `normal_only` status must not be promoted
to calibrated rare-event inference.
