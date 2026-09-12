# Genetically predicted TWAS and PWAS

Source-build feature. These commands test genetically predicted molecular
abundance using GWAS statistics and molecular prediction weights. The existing
numeric-omics association scan remains the workflow for measured expression,
methylation or protein abundance.

## Prepare four aligned inputs

| File | Required columns / layout |
| --- | --- |
| GWAS | `variant ea oa z` |
| Weights | `model variant ea oa weight` |
| Reference | `variant ea oa sd` |
| LD | Square labeled variant correlation matrix, first column containing row IDs |

`weight` is the coefficient per raw effect-allele dosage from a molecular
prediction model. `sd` is the genotype standard deviation from the reference
population. A model may be a gene/tissue, protein, or another genetically
predicted molecular feature. The reference table determines LD orientation.
GWAS z statistics and prediction coefficients are independently aligned to it.
Exact allele swaps reverse signs. No strand complement is guessed; resolve
palindromic/strand-ambiguous variants and genome-build differences upstream.

The model table can contain multiple models. Each model/variant key is unique;
all its variants must exist in reference LD. Every reference variant must have
GWAS statistics. Missing model coverage is an error, so the test does not
silently change the molecular predictor. Use a reference/LD subset appropriate
for the requested locus. Matrix labels are reordered by ID, not by row position.

## Run the example

```powershell
java -jar build/cli/jlinalg-0.3.5.jar twas --gwas examples/xwas/gwas.tsv --weights examples/xwas/weights.tsv --reference examples/xwas/reference.tsv --ld examples/xwas/ld.tsv --joint true --out build/xwas/twas.tsv
java -jar build/cli/jlinalg-0.3.5.jar pwas --gwas examples/xwas/gwas.tsv --weights examples/xwas/weights.tsv --reference examples/xwas/reference.tsv --ld examples/xwas/ld.tsv --out build/xwas/pwas.tsv
```

The shared synthetic example uses two molecular predictors. Its second GWAS
variant is deliberately allele-reversed to exercise alignment. Use your protein
prediction weights for an actual PWAS; the statistical calculation is shared.

For model coefficients w and genotype SDs s, let a = w × s. The test is
`Z = a' z / sqrt(a' R a)`, where R is reference LD. `a' R a` is the variance of
the raw molecular predictor. The first example model has variance 0.05966 and
Z approximately 3.15246 (p = 0.00161903). Scaling every model weight by a positive constant leaves
Z unchanged; reversing every sign reverses Z and preserves its two-sided p-value.

`OUT` reports model ID, weight-row count, Z, unadjusted two-sided p-value, natural log
p-value, and predicted variance. `OUT.metadata.tsv` states weight and LD scale.
The p-value is an association test for the genetically predicted component; it
does not establish a causal molecular effect. Continue with the existing
[colocalization](colocalization.md) and [MR](mendelian-randomization.md) workflows.

## Joint tissue or model test

`--joint true` adds `OUT.joint.tsv`. All models in this invocation form one
prespecified joint family. The command computes their LD-induced correlation
C and the omnibus statistic `Z' C^-1 Z`, with degrees of freedom equal to the
number of linearly independent supplied models. This version requires C to be
positive definite and rejects redundant or numerically ill-conditioned models.
It does not implement S-MultiXcan's rank truncation, model selection or tissue QC.
Split unrelated genes into separate invocations when the scientific question
is a joint test across tissues for one gene.

## Train or import weights

Use [score-train](prediction-scores.md) on a molecular training cohort, with
genotype feature columns, expression/protein as the outcome, and `--alleles`.
The resulting model file can be supplied directly as `--weights`: intercept,
lambda and training metadata columns are ignored by the summary test.
Fit genotype preprocessing inside training folds and validate molecular
prediction outside training. Training expression is quantitative.

External prediction weights can be converted to the same table. Native
PredictDB SQLite and FUSION R-object adapters are not included. Verify that
converted coefficients use raw dosages, the documented allele, and compatible
ancestry/reference LD. Training coefficients on standardized genotypes requires
conversion back to raw dosage units before import.

## Java API and validation

```java
PredictedOmics.Association result = PredictedOmics.test(
    alignedGwasZ, rawDosageWeights, genotypeSd, ldCorrelation);
PredictedOmics.JointAssociation joint = PredictedOmics.joint(
    alignedGwasZ, modelByVariantWeights, genotypeSd, ldCorrelation);
```

Import `org.jlinalg.xwas.PredictedOmics`. LD is row-major. Singular PSD LD is
allowed for a single model if its variance is numerically supported. Invalid
correlations, absent alleles and zero predicted variance fail explicitly.
Dense reference validation is O(variants cubed), so use locus-sized matrices.

Reference: [MetaXcan framework](https://github.com/hakyimlab/MetaXcan).
The [validation report](../xwas-followup-validation.md) distinguishes independent
formula checks from package-wide compatibility claims.
