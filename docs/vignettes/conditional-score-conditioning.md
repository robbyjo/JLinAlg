# Summary-only conditional score analysis

Use this workflow when cohorts have already produced JLinAlg's versioned
conditional-GWAS score exports and the target and conditioning variants share a
complete covariance block. The analysis uses aggregate scores and covariance;
participant-level phenotypes and genotypes remain at each cohort.

This is local, one-step inference at the exported null model. It is useful for
screening a prespecified conditioning set without another data transfer. It is
not a replacement for a cohort-side logistic, Poisson, or Cox likelihood refit
when the conditioning set materially changes a nonlinear null model.

## 1. Plan and export a complete block

The target variants, conditioning variants, and every covariance between them
must be in the same export block. Restricting the genotype input to a locus is
often the clearest way to guarantee that scope. The block size cannot exceed
512 variants.

```powershell
./gradlew.bat executableJar

java -jar build/cli/jlinalg-0.3.5.jar `
  --omics cohort-a-locus.vcf.gz --pheno cohort-a-phenotype.tsv --id IID `
  --formula "case_status ~ age + sex + PC1 + PC2 + <omics>" `
  --family binomial --case-value case --control-value control `
  --conditional-gwas-summary --score-genome-build GRCh38 `
  --score-block-size 128 --out cohort-a-gwas.tsv
```

Repeat the same model contract in each independent cohort. Each successful run
publishes the linked trio consumed below:

- `cohort-a-gwas.tsv`
- `cohort-a-gwas.tsv.score-cov.tsv`
- `cohort-a-gwas.tsv.score-manifest.json`

The JSON manifest is the completion marker. Do not import a partial covariance
file or mix files from different runs.

## 2. Make the cohort manifest

Create a tab-separated file with one row per cohort. Relative paths are
resolved from the directory containing this manifest.

```text
cohort	summary	covariance	manifest
cohort_a	cohort-a-gwas.tsv	cohort-a-gwas.tsv.score-cov.tsv	cohort-a-gwas.tsv.score-manifest.json
cohort_b	cohort-b-gwas.tsv	cohort-b-gwas.tsv.score-cov.tsv	cohort-b-gwas.tsv.score-manifest.json
```

Cohorts are assumed independent. Sample overlap is not modeled, so overlapping
cohorts must not be entered as though their scores were independent.

## 3. Condition and pool

Specify allele-aware keys in the orientation wanted in the output:

```powershell
java -jar build/cli/jlinalg-0.3.5.jar conditional-score `
  --cohorts cohorts.tsv `
  --targets 1:456789:A:G,1:456950:C:T `
  --condition-on 1:455100:G:A `
  --out locus-conditional.tsv
```

For each cohort, JLinAlg performs these operations in order:

1. Match exact forward-strand chromosome, position, and allele pairs. A REF/ALT
   swap is aligned by changing the score sign and the corresponding covariance
   signs only when the exported null is translation invariant: Cox partial
   likelihood or a GLM/OLS null with an explicit `(Intercept)` column. Without
   that guarantee, request the exported orientation; fitted conditioning sets
   across cohorts follow the same rule. Strand complements,
   liftover, and variant normalization are never guessed.
2. Verify the completed version-1 manifest, row-level null-model IDs, fitted
   conditioning-set ID, analyzed sample count, block counts, covariance
   diagonals, and the complete selected block upper triangle.
3. Require compatible genome build, model, family, formula, null covariate
   columns, fitted conditioning set, case/control mapping, analysis-sample and
   genotype/missingness conventions,
   covariance scale, variance method, relatedness adjustment, tie method, and
   tail-calibration declaration across cohorts.
4. Compute the Schur complement separately in each cohort. Singular or
   ill-conditioned conditioning information is rejected.
5. Sum aligned conditional scores and information across the independent
   cohorts and report normal-score inference.

The primary output has one row per requested target with `score_u`,
`score_variance`, the one-step `beta_score`, `se_score`, `z_score`, normal-tail
`p_score_normal`, contributing cohort count, direction string, and the explicit
`local_schur_one_step` scope.

Two linked files are also written:

| File | Purpose |
| --- | --- |
| `locus-conditional.tsv.score-cov.tsv` | Complete upper triangle of the pooled, conditioned target covariance |
| `locus-conditional.tsv.metadata.tsv` | Model contract, allele policy, cohort-independence assumption, and inference boundary |

The ordinary run log defaults to `locus-conditional.tsv.log`; use `--log FILE`
or `--no-log` under the shared CLI logging contract. Existing outputs are never
replaced.

## Missing covariance is not zero

An export block is a coverage unit, not an LD window. If a target is in block 2
and a condition is in block 1, the importer stops with an unknown cross-block
covariance error. It does not insert zero and does not substitute an external LD
matrix. Re-export a locus-sized input with one sufficiently large block, or ask
each cohort to perform a participant-level `--condition-on` refit.

## API use

The same validated file import and conditioning path is available from Java:

```java
import java.nio.file.Path;
import java.util.List;
import org.jlinalg.cli.ConditionalScoreImporter;
import org.jlinalg.settest.ConditionalScoreInference;

var requested = List.of("1:456789:A:G", "1:455100:G:A");
var cohortA = ConditionalScoreImporter.read("cohort_a",
    Path.of("cohort-a-gwas.tsv"),
    Path.of("cohort-a-gwas.tsv.score-cov.tsv"),
    Path.of("cohort-a-gwas.tsv.score-manifest.json"), requested);
var cohortB = ConditionalScoreImporter.read("cohort_b",
    Path.of("cohort-b-gwas.tsv"),
    Path.of("cohort-b-gwas.tsv.score-cov.tsv"),
    Path.of("cohort-b-gwas.tsv.score-manifest.json"), requested);

var fit = ConditionalScoreInference.condition(List.of(cohortA, cohortB),
    List.of("1:456789:A:G"), List.of("1:455100:G:A"));
double conditionalScore = fit.state().scores()[0];
double conditionalInformation = fit.state().information()[0];
```

`ConditionalScoreStudy` is also public for callers that already hold a complete
score vector and covariance matrix. Construction validates dimensions,
allele-aware key uniqueness, positive-semidefinite covariance, and score/null
space consistency.

## Distinguish the three conditioning routes

- General-scan `--condition-on` adds dosages to the cohort null model and refits
  it before export. Use this for likelihood-based nonlinear conditioning.
- `conditional-score` imports the new score schema and uses local Schur-complement
  curvature, cohort by cohort, without participant-level data.
- `mr-estimate --method conditional` remains the existing Gaussian beta/SE plus
  LD approximation. Its input and interpretation are unchanged; it does not
  consume conditional-GWAS score exports.

Focused validation commands:

```powershell
./gradlew.bat test --tests org.jlinalg.settest.ConditionalScoreInferenceTest
./gradlew.bat test --tests org.jlinalg.cli.ConditionalScoreCliTest
```

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Jian Yang et al. (2012) — Conditional and joint multiple-SNP analysis of GWAS summary statistics](../CITATIONS.md#yang-cojo-2012)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
