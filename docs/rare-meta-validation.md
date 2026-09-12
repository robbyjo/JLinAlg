# Rare-variant meta-analysis validation

Validated on 2026-09-12. This report covers the independent-cohort,
quantitative-trait score workflow and its conditional, adaptive, heterogeneous,
diagnostic, and related-sample export extensions described in the
[tutorial](rare-variant-meta-analysis.md). It does not establish parity for every
RAREMETALWORKER model, Raremetal2's exact method, or binary traits.

## Reference and automated checks

The compatibility reference is [statgen/raremetal v4.15.1](https://github.com/statgen/raremetal/tree/v4.15.1),
commit `2c82cfc5710dbd9fd56ef67a7ca5f74772d4e70d`. The checked-in upstream
fixtures and attribution are in
[`src/test/resources/raremetal/PROVENANCE.md`](../src/test/resources/raremetal/PROVENANCE.md).

`RareMetaCliTest` compares single-site, equal-burden, Madsen-Browning weighted
burden, and SKAT results with upstream tutorial outputs. It compares more than
500 informative single sites after allele orientation and excludes malformed
unknown-ALT rows. Tolerances account for the reference's printed precision and
historically rounded covariance files: single p absolute error at most 1e-5,
burden p 2e-5, SKAT p 3e-5. It also exercises indexed rvtests covariance layout,
allele swaps, missing cohorts, minimum cohort counts, exact window boundaries,
unknown alleles in an earlier cohort, and failure on incomplete covariance.

The participant exporter has an independent intercept-only OLS fixture:
RSS=17.5, residual variance=3.5, first-variant score=7/3.5, and cross-covariance
(-4/3)/3.5. BGZF and tabix output are reopened for regional access and all four
meta-analysis tests. Scores retain original phenotype units.

`SummarySetTestsTest` checks analytic single/burden formulas, rank-one SKAT-O,
zero-information burdens, invalid covariance, inconsistent singular scores,
missingness, compensated score summation, and an extreme normal tail whose
ordinary p underflows but log p remains available. It also compares SKAT with
independent R integration to relative tolerance 1e-8 on four rank-three examples.

Regenerate those independent fixtures with:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/rare_meta_reference.R
```

This requires SKAT and CompQuadForm; the recorded versions are R 4.6.1,
SKAT 2.2.5, and CompQuadForm 1.4.4. The generator uses an independent spherical
survival integral for three-dimensional Gaussian quadratic forms. In its extreme
tail case the probability is approximately `1.0640581158717e-22`. Davies returns
zero and Imhof suffers cancellation on this example; neither return is used as
the truth. JLinAlg agrees with the spherical survival integral at the tested
tolerance. This establishes a useful tail case, not universal superiority.

RAREMETAL 4.15.1 does not provide SKAT-O as a selectable test. For SKAT-O,
200,000 correlated Gaussian null simulations are compared with R SKAT's
adjusted result on three moderate-probability cases, within absolute tolerance
0.015. These are different calibration algorithms, so this is a compatibility
check rather than an equality assertion. Default simulation resolution is
`1/(simulations+1)`; its default 10,000 draws cannot resolve genome-wide rare
tails. The explicit analytic option is a moment approximation. The new
`deterministic` option is validated separately below.

Run the complete Java, website, documentation, and executable checks with:

```powershell
.\gradlew.bat check javadoc executableJar --console=plain
```

The original delivery passed its complete check, including website validation,
Javadoc, and executable JAR generation. Expanded validation is described below.
The expanded full run passed **709 tests: 706 passed, three optional native
CHOLMOD skips, zero failures/errors**. Website validation, Javadoc and executable
JAR generation passed. Fresh-output smoke runs also passed both cohort exports,
GRM export, all eight group-test choices with diagnostics, conditional and
deterministic analyses, and both runnable Java examples.

## Advanced methods and rare tails

`SummaryScoreModelsTest` compares conditional blocks with hand-calculated Schur
complements, rejects singular conditioning information and overlapping indices,
and checks a heterogeneous kernel with opposing cohort effects against an exact
chi-square(4) tail. Homogeneous pooling cancels those scores; the heterogeneous
statistic remains 18, as required by the different alternative.

Deterministic SKAT-O is checked against independent **direct polar integration**
in rank two and **chi-square convolution** in rank three. Neither reference uses
the implemented conditional noncentral gamma series or SKAT-O moment matching.

| Reference case | Minimum component p | Adjusted p |
|---|---:|---:|
| Rank 2, moderate | 0.185783307519 | 0.237945171622 |
| Rank 2, rare | 6.33021672719e-13 | 1.18564655285e-12 |
| Rank 2, equal covariance eigenvalues | 2.98872556428e-15 | 5.76925343349e-15 |
| Rank 2, correlated | 0.0586077885131 | 0.0626214859604 |
| Rank 3, moderate | 0.111610225095 | 0.176157013373 |
| Rank 3, rare | 4.18301232668e-21 | 9.01632670699e-21 |

The test requires component relative error <=1e-7 and adjusted relative error
<=3e-6. The implementation targets quadrature error `1e-6 * min-p`, estimates
quadrature error, and bounds positive-series truncation by `1e-10 * min-p`.
It fails explicitly on unresolved spectra, series, or integrals; it does not
replace them with saddlepoint or moment tails. These are algorithmic error
controls, not interval-arithmetic rounding guarantees. The tested cases do not
establish universal convergence, finite-sample calibration, or performance for
large ill-conditioned groups. Supported min-p is >=1e-250 and the noncentral
series is capped at 8,192 terms.

Additional noncentral mixture tests use an exact normal-tail identity for rank
one and independent normal convolution for unequal eigenvalues. The normal
identity was needed to adjudicate cancellation in R's noncentral `pchisq` upper
tail rather than treating R as an oracle. VT is checked against an independently
integrated correlated Gaussian maximum, within five Monte Carlo standard errors;
its output separates selected effects from search-adjusted p-values. A single
threshold reduces exactly to the burden test, and zero-information thresholds
produce an explicit unavailable result.

Regenerate the new base-R-only fixtures with:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/raremetal/advanced-reference.R
```

`RareMetaAdvancedCliTest` exercises all new test names, both leave-out modes,
empty reductions, fixed/random burden pass-through, cohort output, QC column
counts, conditional allele flips on both covariance axes, missing-condition
policies, unknown condition groups, and missing cross-covariance. Serial/no-cache
and two-worker/cached runs must produce byte-identical test and QC tables.

`RareScoreRelatedTest` checks a balanced six-pair random-intercept model against
closed-form ANOVA REML: residual variance `0.86/6`, genetic variance
`3.5 - residual/2`, and independently calculated block inverses for U and V.
Exported scores and covariance must agree to relative/absolute tolerance 2e-6.
It exercises reordered GRM and phenotype IDs, BGZF/tabix round trips, exact HWE,
and fractional-dosage HWE exclusion. The HWE probability comparison uses relative
rather than absolute tolerance, preventing artificial rare-tail p-value floors.

## Dense-block engineering benchmark

Run `gradlew benchmarkRareMetaEngineering --console=plain`. The original synthetic
fixture has four independent cohorts, six groups sharing each dense covariance
block, and 128 or 256 variants per group. Every timed run must produce identical
burden, SKAT, and QC output. It compares one worker with cache disabled, one with
a 16 MiB-per-cohort cache, and two workers with that cache. Results are retained
in a fresh `build/rare-engineering-*/timings.csv` directory. This is an in-process
CLI benchmark with file I/O, repeated numerical validation, and output writing;
it excludes JVM startup and makes no cross-program speed claim.

Median elapsed seconds over three runs on Windows Java 25:

| Variants per group | 1 worker, no cache | 1 worker, 16 MiB cache | 2 workers, 16 MiB cache | Cached worker speed ratio |
|---|---:|---:|---:|---:|
| 128 | 4.6447 | 4.5425 | 2.3731 | 1.91x |
| 256 | 51.7352 | 51.0581 | 26.3398 | 1.94x |

The [recorded timings](../src/benchmark/resources/raremetal/engineering-timings-2026-09-12.csv)
include every repeat. They use plain unindexed synthetic files and repeated
identical regions to exercise cache reuse, not a claim about every BGZF workload.
The larger block remains expensive because eigendecomposition is CPU-bound.

Numerical work dominates these blocks. Cache reuse avoids repeated reads but
does not remove matrix decompositions. Two workers parallelize numerical group
work while shared readers serialize regional I/O. Large heterogeneous kernels
and exhaustive leave-out SKAT-O diagnostics can cost substantially more; the
benchmark does not measure those workloads or imply a genome-scale heap bound.

## Native executable comparison and timing

The pinned C++ reference was built in Ubuntu under WSL using Release CMake,
system zlib, and its pinned cget dependencies. Modern GCC required
`-DCMAKE_CXX_FLAGS="-include cstdint"`; the statistical source was unchanged.
Its obsolete zlib download was replaced by the installed system library. All
recorded reference comparisons disable telemetry with `--noPhoneHome`.

The synthetic benchmark has two independent cohorts of 1,000 participants and
20,000 biallelic sites each, with known score variances and deterministic scores.
Both programs read the same BGZF/tabix inputs. All 20,000 rows agree within the
reference's six-significant-digit output precision:

| Maximum difference | Value |
|---|---:|
| Absolute beta | 4.999411e-7 |
| Absolute SE | 4.315409e-7 |
| Absolute p | 4.999708e-7 |
| Relative p | 4.849571e-6 |

Elapsed seconds for four repeated runs:

| Program | Run 1 | Run 2 | Run 3 | Run 4 |
|---|---:|---:|---:|---:|
| JLinAlg, Windows Java 25 | 0.3461 | 0.3452 | 0.3531 | 0.3372 |
| RAREMETAL, Ubuntu WSL | 0.85 | 0.88 | 0.84 | 0.85 |

The median ratio is about 2.46 in JLinAlg's favor on this workload. Java timing
includes JVM startup; the native timing includes the executable, excluding WSL
launch. RAREMETAL also generates standard plots. Different operating systems,
output precision/features, and a small synthetic single-site workload limit the
comparison. Runs were outside the restricted tool filesystem sandbox, whose
overhead distorted preliminary Java measurements. No speed claim is made for
large groups, SKAT-O, participant export, or genome-scale production data.

To reproduce, build the pinned native executable under
`build/raremetal-reference/build/raremetal`, install `bgzip`/`tabix` in Ubuntu,
then run from the repository root in PowerShell:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/prepare_rare_meta_benchmark.R
.\src\benchmark\prepare-rare-meta-indexes.ps1
.\src\benchmark\benchmark-rare-meta.ps1
```

Use a fresh benchmark directory: compression and indexing refuse existing
outputs. The timing script writes `build/rare-meta-benchmark/timings.csv` and
unique per-run prefixes. Compare a generated Java and native result pair with
`Rscript src/benchmark/r/compare_rare_meta_benchmark.R JAVA.single.tsv NATIVE.meta.singlevar.results`.
Override the scripts' `LinuxWorkspace` parameter if the WSL path differs.

The native reference also consumed JLinAlg's tiny exported example and wrote
matching single-site and equal-burden tables (pooled burden beta approximately
1.58824 and p approximately 0.153007). That two-site/one-group native run exited
with status 1 after writing the tables and did not finalize its log/plots; its
exit issue is unresolved. This is evidence of table interoperability, not a
successful complete native run. The 20,000-site native benchmark completed with
status 0.
