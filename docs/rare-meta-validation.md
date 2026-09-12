# Rare-variant meta-analysis validation

Validated on 2026-09-12. This report covers the initial independent-cohort,
quantitative-trait score workflow described in the
[tutorial](rare-variant-meta-analysis.md). It does not establish parity for every
RAREMETALWORKER model, Raremetal2's exact method, binary traits, conditional
analysis, or heterogeneous-effect meta-analysis.

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
tails. The explicit analytic option is a moment approximation. Deterministic
rare-tail calibration with independent error control remains on the TODO list.

Run the complete Java, website, documentation, and executable checks with:

```powershell
.\gradlew.bat check javadoc executableJar --console=plain
```

The final run passed: 695 tests discovered, 692 passed, three optional native
CHOLMOD tests skipped, and no failures or errors. Website validation, Javadoc,
and executable JAR generation also passed.

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
