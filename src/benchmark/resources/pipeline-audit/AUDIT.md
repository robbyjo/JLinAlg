# Association / GWAS / pipeline / set-test / CLI audit

Audited after dd45e80, 2026-09-08. Direct shared-workspace edits, no commits/push,
no global build/version/site edits. Numerical estimators outside these packages
were inspected as dependencies but not modified. The parent owns final integration.

## Concrete defects and repairs

| Area | Counterexample before repair | Repair / gate |
|---|---|---|
| OLS / GLM score / P3D | Multiplying a valid marker by 1e-8 returned NaN inference | Relative information tolerances; constant-column controls still fail |
| OLS strong signal | y=i+(-1)^i*1e-7, i=0..19: SE 6.16365e-9 vs R 4.07219972e-9 | Direct alternative residual sum of squares |
| REML / adapters | Eight-row intercept null, initial variance .01, maximumIterations=1: nonconverged fit accepted | Refuse nonconverged P3D null; propagate iterative refit status |
| VCF | DS-only row (.25,1.75), no GT: both became NaN | Read DS before testing GT availability; shuffled sample order asserted |
| QC | Dosages 1-1e-8,1,1+1e-8: variance zero and monomorphic exclusion | Centered variance; positive variance is not monomorphic; infinity rejected |
| Alignment / transforms | Duplicate/negative row indices accepted; shortened transformed vector zero-filled trailing sample | Validate index uniqueness and exact transform length; defensive input copy |
| Accounting | signal/filtered/collinear produced filtered,ok,ok,failed events | One source-ordered outcome per item; failures not counted as successful tests |
| Annotations | Absent key with selected annotation columns threw NPE | Empty strings for unmatched cells |
| BH records | Embedded newline caused EOF while assigning q-values | Logical quoted records, including quoted headers |
| BH publication / resources | Non-overwrite atomic move could clobber concurrent output; unbounded open sort runs | Non-replacing move, 64-way multipass merge, cleanup and scratch retention rules |
| CLI paths / resume | --out input --overwrite destroyed input; --resume accepted arbitrary existing file | Preflight aliases including logs/manifest/partial; reject unverifiable resume |
| Standalone CLI outputs | Regression/MR estimator output was silently truncated without an overwrite option | Fresh-output preflight and CREATE_NEW for result tables |
| Set-test scale | Weight 1e-8 discarded the positive kernel eigenvalue | Relative spectrum checks, scale-normalized tails, relative burden checks |
| Rank-two mixture | lambda=(1,10), Q=20: fallback p=.17618236 vs R .16908111409965279 | Nonoscillatory angular integral; critical values invert component tails |
| Independent review of new integral | lambda=(1,1e-12): q50/q200 had 25%/38.5% tail errors from absolute tolerance | Factor out exp(-Q/(2*lambdaMax)); extreme-tail/scale tests added |
| Pasteur rank-three review | lambda=(1,.5,.25), q100: Imhof returned 1.207e-10 although chi-square(3) upper bound is 1.554e-21 | Positive gamma-mixture series gives R-matching 2.5052827321161243e-23; explicit relative truncation bound and resource failure |
| Pasteur extreme predictor units | x*1e-160 gave infinite SE/NaN p; x*1e160 falsely nonestimable | Normalize before imputation/projection; map beta and SE back; unrepresentable direct OLS exports recorded as failures |
| Pasteur padded-rank quantile | Ten eigenvalue slots with one positive value, p=nextDown(1): q4.77e-17 vs exact1.936e-32 | Positive-rank exact branches, relative bracketing, direct lower-tail inversion near p=1 |

The independent review probe now gives tail/chi-square(1) ratios
1.0000000000046225 at Q=50 and 1.000000000000637 at Q=200.
The near-rank-one gates also cover Q=700, with reference p=2.9902269751246203e-154,
and common eigenvalue/statistic scales 1e-100, 1, and 1e100.

## Validation and measured performance

After the final Pasteur repairs, the complete isolated owned-package run passed
83/83 tests (175.469 seconds), including the new extreme-scale, bounded-series,
and near-one quantile regressions. Pasteur independently recompiled the three
changed production files and cleared all three original findings, plus additional
near-one ellipsoid bounds, p=1e-100 inversion, scales 1e+/-200, and the BH/alignment
probes. His historical report and final independent evidence are in
`build/audit/GAUSS_CROSS_REVIEW.md`. Both OLS/BH benchmark runners assert frozen R checksums;
the BH runner creates its scratch directory on a clean checkout.
Logs are in build/pipeline-audit/all-owned-tests.txt and
build/pipeline-audit/final-output-kernel-tests.txt. Parent integrated tests are
separate evidence and must be rerun after the final source snapshot.

R-generated fixtures are in src/test/resources/pipeline-audit. The 128-marker
400-sample full OLS scan includes deterministic missing values, mean imputation,
intercept/covariate adjustment, alternative RSS, residual DF and t probabilities.
Max absolute discrepancies after predictor normalization: beta 2.49801e-16,
SE 1.11022e-16, p 1.66534e-15.
BH fixtures explicitly exclude failed/NaN rows from R p.adjust's input family.
Kernel references use base R integrate over the independent polar identity and
a separately generated positive gamma-mixture fixture with mass/mean checks.

Raw CSVs alongside this report retain negative-repetition warmups and every timed
repetition, with consumed checksums. Baseline classes were compiled from exact
`git show dd45e80:src/main/java/org/jlinalg/{association/FastOlsAssociation.java,cli/ExternalBh.java}`
sources into build/pipeline-audit/baseline-classes; all other dependencies match the
current benchmark environment. This is important because the parent's build updates
build/classes/java/main during integration; that directory alone is not a baseline.

| Workload | Baseline Java median ms | Repaired Java median ms | R median ms |
|---|---:|---:|---:|
| OLS 400 samples x 128 markers | .8204 | .7514 | 5.2000 |
| Disk BH, 20,000 rows | 91.1752 | 78.1121 | 1.0000 (memory only) |
| Rank-two tail Q=20, lambda=(1,10), per call | not accuracy-equivalent | .0293035 | .022000 |
| Rank-three tail Q=100, lambda=(1,.5,.25), per call | not accuracy-equivalent | .114658 | 3.2000 |

OLS: 30 Java warmups then 9 fits; R 5 warmup batches then 9 timed batches of 25
scans to resolve the Windows clock. Both fit the same alternative OLS estimand;
Java caches the covariate projection whereas R uses a fresh QR per marker. OLS
checksums agree at 8075.53125328881. Java is about 6.9x faster than this R loop
in this final run. The retained baseline timing is historical, not a simultaneous
paired comparison: do not attribute its difference solely to the last accuracy fix.

BH: 3 warmups, 7 timed runs, whole write/sort/adjust/publish operation. R uses
20 in-memory adjustments per timing batch and excludes I/O. Checksums agree at
176709.280009885; current Java is not faster than R's memory-only work.

Rank two: 5 warmup batches, 9 timed batches, Java 200 and R 5000 integrations per
batch. Checksums are respectively 33.81622281992951 and 845.405570498212
(different repetition counts, same p). Java is slower than R on this fixture.
No optimizer nonconvergence contributes to these timing claims.

Rank three: 5 warmup batches, 9 timed batches, Java 100 and R 25 tail evaluations
per batch. Both use the same positive gamma-mixture identity, and both generate
coefficients inside every timed call. Java adaptively stops at its relative PGF
remainder bound; R fixes 600 terms and independently asserts the same 1e-12 bound
for this q=100 workload. Java/R per-batch checksums are respectively
2.5052827321161185e-21 and 6.26320683029031e-22 (different counts, same tail).
Java was about 28x faster on this bounded fixture, not a guarantee across spectra.
Raw final Java/R timings were refreshed together after the 83-test run; all
reference-checksum gates passed. The registered benchmark task includes both ranks.

## Reproduction (PowerShell, repository root)

```powershell
.\gradlew test --tests 'org.jlinalg.association.*' --tests 'org.jlinalg.gwas.*' --tests 'org.jlinalg.pipeline.*' --tests 'org.jlinalg.settest.*' --tests 'org.jlinalg.cli.*'
.\gradlew benchmarkPipelineOlsAudit benchmarkPipelineBhAudit benchmarkPipelineKernelAudit --no-parallel

# Optional: regenerate the independent base-R fixtures and warm timing evidence.
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/pipeline-audit/generate-reference.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/pipeline-audit/generate-rank-three-reference.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/benchmark/r/pipeline_audit_benchmark.R
```

The ordinary Gradle test report is `build/reports/tests/test/index.html`.
No ignored local launcher or copied JUnit dependency directory is required.
Each benchmark task prints CSV warmups/timings and asserts its frozen reference
checksum; the kernel task includes both rank-two and rank-three workloads.
Run timing tasks sequentially on an otherwise idle host. During shared-agent
editing the parent alone coordinates Gradle; the above commands are the normal
post-integration/clean-checkout reproduction interface.

### Executable CLI smoke after parent builds the final executableJar

Use a fresh smoke directory (the --overwrite here intentionally reruns only the
named build output, never the source fixtures):

```powershell
java -jar build/cli/jlinalg-0.3.0.jar --help
java -jar build/cli/jlinalg-0.3.0.jar mr-estimate --help
java -jar build/cli/jlinalg-0.3.0.jar beta-regression --help
java -jar build/cli/jlinalg-0.3.0.jar --pheno src/test/resources/pipeline-audit/smoke-pheno.tsv --omics src/test/resources/pipeline-audit/smoke-omics.tsv --id id --formula 'y ~ <omics>' --model ols --backend cpu --threads 1 --omics-type generic --annot src/test/resources/pipeline-audit/smoke-annotation.tsv --annot-cols all --out build/pipeline-audit/smoke.tsv --overwrite
```

Expected: exit 0; header plus three rows in order signal/constant/noise, statuses
ok/failed/ok, missing annotations blank, two finite BH values, manifest counts
source=3/tested=2/failed=1. The phenotype fixture deliberately shuffles sample IDs.
Repeat without --overwrite: nonzero exit and unchanged outputs. Passing the
phenotype path as --out or --log is rejected before writes, even with --overwrite.
Passing --resume with an existing arbitrary result is nonzero, not a false success.
Run the package-level VCF test for DS-only data and requested-order validation.

## Limits and deliberately unchanged semantics

All owned package sources and three benchmark runners compile in isolation.
Each runner explicitly rejects nonfinite checksums (NaN must not bypass a
tolerance comparison). Final raw timing evidence is versioned alongside this
report; ignored build-local files are not required for normal reproduction.

- Higher-rank positive mixtures use nonnegative gamma-series terms, not
  oscillatory Imhof cancellation or a Satterthwaite fallback. For coefficient
  PGF G(z) and 1<z<1/max(r_i), omitted mass after K is bounded by G(z)/z^(K+1).
  That bound must be at most 1e-12 of the computed probability (with a tighter
  gamma-CDF factor for lower tails). The limit is 8192 terms, O(K^2+rank*K)
  work and O(K+rank) storage. Ill-conditioned/unrepresentable spectra or failure
  to satisfy the bound throw an explicit numerical exception. The bound controls
  truncation, not interval-arithmetic roundoff; not every positive spectrum is
  guaranteed to finish. Analytic SKAT-O remains moment-matched, independently of
  the component mixture-tail algorithm; no exact finite-sample SKAT-O claim.
- Mixture probabilities below representable range retain the existing positive
  Double.MIN_VALUE floor; no new log-tail API is introduced here.
- P3D is fixed-null variance-component inference; GLM and GAM prepared scans
  remain score approximations. Raw matrix APIs cannot verify biological IDs
  beyond their explicit index/ID contracts. No estimator internals outside
  ownership were changed.
- BH excludes failed/invalid p-values from its denominator. A user prespecifying
  all attempted hypotheses must deliberately decide how failures enter that family.
- Resume/checkpoint authentication is not implemented. A partial file is retained
  for recovery, not automatically trusted. Main CLI path checks are not a security
  boundary against an adversary concurrently changing filesystem links.
- LD archive extraction was inspected read-only: expected basenames are mapped to
  configured targets, duplicate entries rejected, nonregular entries skipped.
  No live downloads/network security audit or large hostile-archive fuzzing was run.
- No performance generalization across backends, cohorts, storage, or rare tails.

## Changed paths

Sources: association/{AssociationModels,FastGlmAssociation,FastOlsAssociation,
ParallelAssociationEngine}; gwas/RemlAssociationScanner;
pipeline/{SampleAlignment,StreamingAssociationPipeline,StreamingOmicsAssociationPipeline,
VariantFilters,VariantStatistics,VcfVariantSource}; settest/{QuadraticFormDistribution,
RemlSetTestNullModel,SetTests}; cli/{AnnotationLookup,CliOptions,ExternalBh,JLinAlgCli,
PipelinePaths,LdClumpCli,MrEstimatorCli,MrInstrumentCli,MrXwasCli,RegressionCli}.
All are Java files under src/main/java/org/jlinalg.

Five audit test classes under matching src/test/java packages; all files under
src/test/resources/pipeline-audit; the three Pipeline*AuditBenchmark.java classes;
src/benchmark/r/pipeline_audit_benchmark.R; this report and raw timing CSVs;
docs/vignettes/association-gwas-twas.md and docs/vignettes/xwas-mr-pipeline.md.
