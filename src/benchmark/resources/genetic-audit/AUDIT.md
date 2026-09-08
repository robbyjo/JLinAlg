# Post-dd45e80 coloc, SuSiE, genetics, and MR audit

Date: 2026-09-08. Parent owns integration/release. Implementation used isolated
compilation during concurrent edits; the clean-checkout reproduction below uses
the normal Gradle workflow. No commits/pushes, validated conditional-path
rewrites, or global build/release/site edits were performed by this auditor.

## Reproduce

From the repository root in PowerShell:

```powershell
.\gradlew.bat check
.\gradlew.bat benchmarkGeneticAudit
```

Gradle fetches the declared dependencies, compiles the normal source sets,
and supplies the regression/benchmark classpaths. No pre-existing `build/`
classes, manually copied JUnit jars, or ignored dependency directories are
required. The committed R-generated fixture is sufficient for these commands;
R is not required to run the regression tests or Java accuracy/timing gates.
Run timings sequentially, without concurrent builds or other benchmarks.

To regenerate the independent base-R fixture and capture fresh R timings:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-genetic-audit-reference.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/r-reference/generate-genetic-final-fixes.R
.\gradlew.bat check
.\gradlew.bat benchmarkGeneticAudit
```

The R script prefixes `.libPaths` and needs no additional R packages. A
PowerShell convenience wrapper runs the same Gradle commands:
`& src/test/resources/r-reference/run-genetic-audit.ps1`. Add `-RegenerateR`
to regenerate first; override its R executable with `-Rscript 'path/to/Rscript'`.
The wrapper does not install R or packages. During shared-workspace integration,
the parent coordinates execution of the full Gradle commands.

Sources: `src/test/resources/r-reference/generate-genetic-audit-reference.R`,
`genetic-audit-reference.properties` beside it, and
`src/benchmark/java/org/jlinalg/genetics/GeneticAuditBenchmark.java`.
[Final-fix measurements](final-fix-evidence.txt) include checksums and call counts;
[earlier measurements](raw-warm.txt) are retained as historical evidence. Benchmarks
consume the sum of index-weighted squared outputs, avoiding zero checksums
from centered GRMs. Different repetition counts mean compare checksum/calls,
not total checksums. Every Java single-effect fit must converge. The recorded
raw measurements came from the earlier isolated development runner, not a
newly measured clean-checkout Gradle run; that runner is no longer required.

## Reproduced defects and repairs

| Area | Counterexample before repair | Repair/reference gate |
| --- | --- | --- |
| Coloc disjoint support | BF rows `[100,-Inf]`, `[-Inf,100]` gave H0=1 | H3=1 to machine precision; enumerate distinct configurations in R |
| Coloc cancellation | Both rows `[1000,950]` gave H3=0 | Direct O(P) ordered-distinct log accumulation retains H3/H4=`0.004 exp(-50)/(1+exp(-100))` |
| Coloc weighted SNP posterior | H4 hypotheses used weights but conditional SNP probabilities did not | Product weights appear in both; R enumeration and weight-scale/permutation checks |
| Coloc weight range | Finite weights near 1e308 overflowed their normalizer | Log normalization, invariant to positive vector rescaling |
| SuSiE summary overflow/PIP | z=`1e200` gave PIPs `[.5,.5]` and zero effects | `hypot` transformation gives saturated beta .9519230769 and weak PIP 3.435e-21; log-complement PIP accumulation |
| SuSiE iteration state | One sweep returned beta 1/3 at sigma²=1 but exposed sigma²=.7676768 | Exhaustion preserves the variance used by the returned posterior/objective, with convergence false |
| SuSiE validation/convergence | Indefinite LD accepted as converged; large negative ELBO step met old stopping rule | External PSD/symmetry checks, finite/integer scales and unique names, joint sufficient-statistic coherence, absolute objective-change criterion |
| GRM | Dosages `[0,1,2]` rejected the centered middle sample's zero diagonal | Allow zero rows/singular PSD; check external covariance with scaled pivoted Cholesky |
| IVW/Egger scaling | IVW rejected common 1e-200 scaling; exact `y=2x+3` near x=1e8 returned slope 4 | Scaled accumulation and anchored weighted centering; analytic slope/SE and unit-invariance gates |
| Steiger | Squaring finite beta/SE at extreme units produced NaN R² | Stable standardized t/hypot conversion; matched estimand across unit scalings |
| Winner's curse | t=40.01, c=40 returned 2.96074 after selection probability floor | Exact conditional score root .107180426634; independent R Brent root .107180426636 |
| MR overdispersion | Eight symmetric residuals with second moment 1 and known variance 1 returned tau²=.5 | Truncate the averaged excess, not each excess; tau²=0 within rounding; beta sandwich uses actual adjusted-score contributions |
| Mixture likelihood | Outcomes `[-400,0,400]`, x=1, sy=1 gave logLik=-1382.4702 | Log-sum-exp likelihood -158418.715662496; same grid/profile curvature as base R |
| Multivariable strength | `conditionalFStatistics()` was mean marginal z² | New honest `marginalFStatistics()`; deprecated old accessor explicitly throws instead of reporting fictitious conditional strength |
| Unphased clumping LD | Four corner homozygotes plus 100 double heterozygotes gave r²=0; subsequent multi-start EM could time out or stop at a nonmaximum | All admissible cubic stationary phases and endpoints, compared by observed likelihood; analytic and independent R polyroot gates |

All 64 tests passed in the final isolated four-package run. This includes the existing
574×1,001 N3 susieR fixture at unchanged PIP/coefficient 2e-10 and objective
2e-8 tolerances, existing coloc package fixtures, and prior conditional,
generalized/overlap-aware MR, allele alignment, confidence and plot tests.
There are new dedicated boundary tests in every owned package and a base-R
reference test covering all 69 coloc, 96 single-effect, and nine GRM outputs.

## Warm timing and accuracy

Five batches after warming, same fixed inputs, sequential R then Java on the
shared Windows host (Java 25, 32 available processors; R 4.6.1):

| Workload | Java ms/call | Base R ms/call | Accuracy gate |
| --- | ---: | ---: | --- |
| 64-SNP coloc, one pair | 0.01244 | 0.2 | max scaled absolute error 5.55e-17 |
| 32-SNP one-effect posterior | 0.02307 | 0.006 | 8.88e-16 |
| Selected-normal MLE t=40.01, c=40 | 0.02237 | 0.12 | 2.17e-12 across four MLE cases |
| Three-variant/three-sample GRM | 0.000813 | 0.0035 | zero at the recorded precision |
| 41-founder LD stationary phase | 0.003293 | 0.030 | relative r² error 4.22e-15 |

These are the refreshed paired measurements after the final SuSiE conjugate
arithmetic repair, not the earlier implementation snapshot.
The LD row is a subsequent paired additive run on the same host: Java includes
reflection dispatch, genotype counting, stationary roots and likelihood selection;
R starts from the same pre-counted 3×3 table and evaluates all `polyroot` candidates
and endpoints. Neither timing includes BED I/O. This is not a PLINK speed comparison.

These are **base-R equation/reference timings, not fresh package benchmarks**.
R coloc explicitly enumerates the P×P configuration matrix; Java accumulates
the equivalent sums in O(P). R's one-effect reference is closed form; Java
validates LD and runs two IBSS sweeps to convergence, so it is slower here.
R GRM times the pre-standardized crossproduct, whereas Java times dosage
filtering, imputation, standardization and construction too. These tiny cases
do not establish large-region throughput or a general GRM/SuSiE speedup.
`coloc` and `susieR` were not loadable in this session; their already-bundled
reference fixtures were retained and tested, with no new package-speed claim.

## Actual limits and integration notes

- New external LD/GRM PSD validation costs O(P³) in the full-rank worst case;
  known individual-level crossproducts skip that check. Singular PSD is allowed,
  but the check is numerical (1e-10 scaled pivot tolerance), not a PSD repair.
- External ancestry/allele alignment cannot be inferred from a numeric LD
  matrix. Joint coherence is enforced for claimed exact sufficient statistics,
  not for approximate external-LD RSS summaries.
- Coloc overlap trimming and per-trait leads remain uniform-BF-based. Weights
  alter hypothesis and H4-conditional posterior evidence, not the upstream fit.
- `conditionalFStatistics()` is an intentional compatibility correction and
  should be called out in release notes. No true multivariable conditional F
  estimator was invented without a cross-exposure covariance model.
- Multivariable Egger uses caller-chosen joint allele orientation; it does not
  silently impose an exposure-increasing reference. Generalized/overlap-aware
  estimators and validated conditional association methods were not rewritten.
- Huber adjusted-profile MR retains a nonrobust moment dispersion estimate and
  fixed-dispersion beta sandwich, not full jointly calibrated R mr.raps parity.
  The PRESSO-style method remains an analytic diagnostic, not R simulation
  calibration. Neither is advertised as a complete package-equivalent method.
- Mixture invalid variance remains fixed at 0.01 in outcome units and valid
  mixing weights are grid-restricted to .05–.95. Curvature inference is grid
  local, not an exact profile CI; unresolved curvature now fails explicitly.
- Haplotype phase fitting evaluates all real stationary-root candidates and
  endpoints rather than iterating EM. Multiple maximizing phases may remain
  unidentified. The existing one-variant ieugwasr pass-through convention stays.
- Winner correction requires a selected observation and threshold <=1e6. It
  corrects an exposure effect, not a conditional instrument-association p-value.
- None of these tests prove universal absence of numerical defects. Raw MR
  double arithmetic can still exceed representable range outside the tested
  scalings; the conditional and statistical assumptions remain essential.

## Final independent-cross-review closure

Gauss's three additional original findings were reproduced, repaired, and
independently rechecked after fresh compilation. All are closed:

- **LD global phase:** the ordinary 41-founder table
  `[[4,0,2],[3,21,4],[3,3,1]]` previously aborted at one EM timeout. It now
  returns r²=0.04048732776637871 and the public BED clump excludes the second
  SNP at 0.01. The table `[[2,2,1],[1,13,3],[1,3,0]]` has a stationary but
  nonmaximizing midpoint; its global maxima give r²=1/675, so the second SNP
  is correctly excluded at the default 0.001 threshold. Permanent tests run
  both through the actual BED/public clump entry point.
- **SuSiE tiny prior:** scalar `X'X=1e-200`, `X'y=1e-100`, `y'y=1`, `n=2`,
  prior variance `1e-200` gives mean `1e-300`, with convergence true. The
  smaller-variance conjugate branch avoids discarding the prior in an
  underflowed shrinkage ratio; exponent-separated products retain the mean.
  Both effect signs are tested.
- **IVW finite reconstruction:** exposure `[1e-200,1e-200]`, outcome
  `[1e200,1]`, outcome SE `[1e200,1]` now gives beta and SE `1e200`, rather
  than an exception from an overflowing intermediate outcome/exposure scale
  ratio. Reconstruction separates binary exponents.

The LD algorithm forms the cubic `F(t)=t(A(t)+B(t))-A(t)` for coupling and
repulsion products. Its derivative roots partition `[0,1]` into monotone
intervals; sign-changing roots are bisected, repeated-root candidates and
endpoints are retained, and observed genotype likelihood selects the global
candidate. It neither skips failed EM starts nor assumes every stationary
point is a maximum. Independent base R uses `polyroot` plus full observed
likelihood evaluation, not Java's root-partition algorithm. The original
bounded 2,000-table reviewer probe completed with no failures.

Permanent regeneration and fixture:
`src/test/resources/r-reference/generate-genetic-final-fixes.R` and
`genetic-final-fixes.properties` beside it. The normal Gradle check covers
the four new regressions; the wrapper's `-RegenerateR` refreshes both reference
scripts. [Final raw evidence](final-fix-evidence.txt) records accuracy, tests,
independent closure, and the refreshed warm timings. Production is frozen for
parent integration; no remaining finding from this final review is open.

## Changed paths for integration

Production: `coloc/ColocSusie.java`, `coloc/ColocSusieResult.java`;
`susie/Susie.java`, `susie/SusieOptions.java`;
`genetics/GeneticCovarianceValidation.java`, `GenomicRelationshipMatrix.java`,
`PlinkBedLdClumper.java`; `mr/MendelianRandomization.java`,
`WinnerCurseCorrection.java`, `SteigerFiltering.java`,
`MultivariableMendelianRandomization.java`, `MultivariableMrResult.java`,
`RobustMendelianRandomization.java`, `ContaminationMixture.java`, all under
`src/main/java/org/jlinalg/`.

Tests: `ColocSusieTest`, new `ColocAuditBoundaryTest`, new
`SusieAuditBoundaryTest`, `PlinkBedLdClumperTest`, new
`GeneticCovarianceAuditTest`, new `GeneticAuditRReferenceTest`, new
`MrAuditBoundaryTest`, in their corresponding `src/test/java/org/jlinalg/`
packages. Reference/benchmark paths are listed above and beside this report.
Source vignettes: `docs/vignettes/colocalization.md`, `susie-and-sem.md`
(SuSiE section only), and `mendelian-randomization.md`.
