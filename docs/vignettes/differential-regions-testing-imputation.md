# Differential, regional, multiple-testing, and imputation workflows

These source-build workflows cover four common points after an omics matrix has
been quality controlled: differential analysis, spatial EWAS aggregation,
prespecified multiplicity correction, and uncertainty-aware missing-data
inference. They are deliberately separate commands because count likelihoods,
region definitions, testing families, and imputation uncertainty are different
statistical contracts.

Build the executable before using these examples:

```powershell
.\gradlew.bat clean check executableJar --no-daemon --no-parallel
$jar = 'build/cli/jlinalg-0.3.6.jar'
```

The commands below are included in v0.3.6 and later.

## Empirical-Bayes differential analysis

The feature matrix is feature-by-sample with a row identifier in column one.
The phenotype table supplies the sample ID, a two-level group, and optional
numeric covariates. Samples are intersected in matrix order; the first aligned
group level is the reference and the second is tested.

```text
feature_id  s1  s2  s3  s4  s5  s6
gene_a      10  12  11  30  33  29
gene_b      50  51  49  52  48  50
```

```text
sample_id  condition  age
s1         control    40
s2         control    41
s3         control    39
s4         case       40
s5         case       41
s6         case       39
```

Continuous intensities use moderated Gaussian linear models:

```powershell
java -jar $jar differential --method limma `
  --omics methylation.tsv --pheno samples.tsv --id sample_id `
  --group condition --covariates age --out build/differential-limma.tsv
```

Counts have two distinct choices. `voom` calculates normalized log counts,
estimates a mean/variance precision trend, fits weighted linear models, and
moderates their residual variances. `negative-binomial` uses median-ratio size
factors, feature dispersions shrunk on the log scale, and log-link NB IRLS.

```powershell
java -jar $jar differential --method voom `
  --omics counts.tsv --pheno samples.tsv --id sample_id `
  --group condition --out build/differential-voom.tsv

java -jar $jar differential --method negative-binomial `
  --omics counts.tsv --pheno samples.tsv --id sample_id `
  --group condition --out build/differential-nb.tsv
```

The primary table contains the effect, log2 fold change where applicable,
standard error, statistic, degrees of freedom, raw and moderated variance or
dispersion, p-value, BH FDR, convergence, and iterations. Companion files
record every sample's size factor and the fitted prior. The BH family contains
all fitted features, not only rows that passed a significance filter.

These compact estimators match the documented contracts and checked fixtures;
they are not full reimplementations of every limma, edgeR, or DESeq2 option.
There is no automatic contrast search, outlier replacement, independent
filtering, quasi-likelihood empirical Bayes layer, or fold-change prior.

## Region-level EWAS

Region input is one row per probe with explicit build and coordinate columns:

```text
probe_id  chromosome  position  effect  p
cg001     chr1        101230    0.08    0.002
cg002     chr1        101370    0.05    0.010
cg003     chr1        101540    0.07    0.004
```

```powershell
java -jar $jar ewas-regions --input probe-results.tsv `
  --genome-build GRCh38 --max-gap 500 --min-probes 3 `
  --correlation-length 200 --out build/ewas-regions.tsv
```

Probes are canonicalized and sorted by chromosome and position. Duplicate
coordinates are rejected. A region is a nonoverlapping, same-chromosome run in
which adjacent probes are no farther apart than `--max-gap`. Signed two-sided
probe Z statistics are summed with covariance
`R(i,j) = exp(-abs(position_i-position_j)/correlation_length)`. This makes the
spatial-dependence assumption visible instead of treating neighboring probes as
independent. Region BH is calculated over every retained region.

Outputs retain the genome build, start/end, complete probe list, probe count,
probes per kilobase, maximum observed gap, correlation length, mean effect,
spatial Z, p-value, and FDR. The command does not lift coordinates or infer a
build. Change-direction mixtures can cancel by design; use a separately
prespecified unsigned regional hypothesis if cancellation is not the estimand.

## Adaptive and hierarchical multiple testing

IHW-style cross-weighting needs a p-value, an informative covariate, and an
independent fold. The covariate must be independent of the p-value under every
true null. The CLI makes the analyst acknowledge that design condition:

```powershell
java -jar $jar multiple-test --method ihw --input tests.tsv `
  --id id --p p --covariate mean_normalized_count --fold fold `
  --null-independent true --bins 5 --out build/tests-ihw.tsv
```

Weights for one fold are learned only from other folds, normalized to mean one
within the held-out fold, and passed to weighted BH. This is a compact
cross-fitted bin-weighting estimator, not the full convex IHW optimizer.
FDR interpretation requires independent folds and the null-independence or
appropriate positive-dependence assumptions documented for weighted BH.

For a prespecified forest, every internal and leaf hypothesis is supplied:

```text
id                 parent       p       status
gene_G                           0.001   ok
gene_G:tissue_A    gene_G       0.010   ok
gene_G:A:trait_1   gene_G:tissue_A  NA  failed
```

```powershell
java -jar $jar multiple-test --method hierarchy --input hierarchy.tsv `
  --id id --parent parent --p p --status status --alpha 0.05 `
  --out build/hierarchy-adjusted.tsv
```

Leaf-count weights are normalized across every tested node, then weighted
Bonferroni p-values are gated through all ancestors. This controls familywise
error under arbitrary dependence and ensures that a child cannot be rejected
when its parent is not. It is intentionally not labeled hierarchical FDR.
Failed rows remain in the complete prespecified family with p=1 and appear in
the output and metadata.

## Multiple imputation and pooled inference

Specify every modeled variable and its compatible imputer. Missing tokens are
blank, `NA`, `NaN`, or `.`.

```powershell
java -jar $jar multiple-impute --input incomplete.tsv --id sample_id `
  --types age:continuous,case:binary,site:categorical `
  --imputations 20 --iterations 20 --donors 5 --seed 1709 `
  --out build/imputed/cohort
```

Every conditional update first bootstraps its observed rows, propagating
model-parameter uncertainty. Continuous variables use predictive mean matching
within that bootstrap donor sample, randomizing ties for each missing cell.
Binary and categorical responses use one joint multinomial logistic model;
nominal predictor categories use dummy columns. This bootstrap MI method uses
weak ridge regularization for separated bootstrap samples and stops if a
conditional fit fails to converge. Each completed data set uses an
independent deterministic stream derived from the recorded seed. The command
writes `OUT.impN.tsv`, chain mean/variance diagnostics for originally missing
cells, and run metadata. Observed cells are never replaced.

Fit the same prespecified analysis separately to every completed data set and
assemble one row per parameter and imputation:

```text
parameter  imputation  estimate  variance
case       1           0.42      0.010
case       2           0.47      0.012
case       3           0.39      0.011
```

```powershell
java -jar $jar mi-pool --input estimates.tsv --complete-df 120 `
  --out build/imputed/pooled.tsv
```

Pooling reports the mean estimate, within- and between-imputation variance,
total variance, standard error, relative variance increase, fraction of missing
information, statistic, p-value, and Barnard-Rubin finite-sample degrees of
freedom. Omit `--complete-df` only for a justified large-sample analysis.
Multiple imputation assumes an adequate imputation model and missing at random
conditional on its predictors; reproducible draws do not make that assumption
true. Deterministic mean imputation is not used and is not equivalent.

## Java API map

The command implementations are thin file adapters over these public APIs:

```java
DifferentialFit fit = EmpiricalBayesDifferential.fitVoom(
    countsByFeature, designBySample, contrast);

List<EwasRegion> regions = RegionLevelEwas.scan(
    probes, "GRCh38", RegionLevelEwas.Options.defaults());

IndependentHypothesisWeighting.Result weighted =
    IndependentHypothesisWeighting.adjust(p, covariate, folds,
        IndependentHypothesisWeighting.Options.defaults());

MiceResult imputed = MiceImputer.impute(data, variableTypes,
    new MiceOptions(20, 20, 5, 1709L, 1e-6));
RubinPooling.Estimate pooled = RubinPooling.pool(beta, variance, 120);
```

See the [validation report](../inference-workflows-validation.md) for frozen
Bioconductor versions, exact numerical gates, reproducibility commands, and
the supported-scope audit.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Gordon K. Smyth (2004) — Linear models and empirical Bayes methods for assessing differential expression in microarray experiments](../CITATIONS.md#smyth-limma-2004) — [PMID: 16646809](https://pubmed.ncbi.nlm.nih.gov/16646809/)
- [Charity W. Law, Yunshun Chen, Wei Shi, and Gordon K. Smyth (2014) — voom: precision weights unlock linear model analysis tools for RNA-seq read counts](../CITATIONS.md#law-voom-2014) — [PMID: 24485249](https://pubmed.ncbi.nlm.nih.gov/24485249/) · [PMCID: PMC4053721](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053721/)
- [Mark D. Robinson, Davis J. McCarthy, and Gordon K. Smyth (2010) — edgeR: a Bioconductor package for differential expression analysis of digital gene expression data](../CITATIONS.md#robinson-edger-2010) — [PMID: 19910308](https://pubmed.ncbi.nlm.nih.gov/19910308/) · [PMCID: PMC2796818](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/)
- [Michael I. Love, Wolfgang Huber, and Simon Anders (2014) — Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](../CITATIONS.md#love-deseq2-2014) — [PMID: 25516281](https://pubmed.ncbi.nlm.nih.gov/25516281/) · [PMCID: PMC4302049](https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/)
- [Brent S. Pedersen, David A. Schwartz, Ivana V. Yang, and Katerina J. Kechris (2012) — Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values](../CITATIONS.md#pedersen-combp-2012) — [PMID: 22954632](https://pubmed.ncbi.nlm.nih.gov/22954632/) · [PMCID: PMC3496335](https://pmc.ncbi.nlm.nih.gov/articles/PMC3496335/)
- [Timothy J. Peters et al. (2015) — De novo identification of differentially methylated regions in the human genome](../CITATIONS.md#peters-dmrcate-2015) — [PMID: 25972926](https://pubmed.ncbi.nlm.nih.gov/25972926/) · [PMCID: PMC4429355](https://pmc.ncbi.nlm.nih.gov/articles/PMC4429355/)
- [Yoav Benjamini and Yosef Hochberg (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](../CITATIONS.md#benjamini-hochberg-1995)
- [Nikolaos Ignatiadis, Bernd Klaus, Judith B. Zaugg, and Wolfgang Huber (2016) — Data-driven hypothesis weighting increases detection power in genome-scale multiple testing](../CITATIONS.md#ignatiadis-ihw-2016) — [PMID: 27240256](https://pubmed.ncbi.nlm.nih.gov/27240256/) · [PMCID: PMC4930141](https://pmc.ncbi.nlm.nih.gov/articles/PMC4930141/)
- [Nicolai Meinshausen (2008) — Hierarchical testing of variable importance](../CITATIONS.md#meinshausen-hierarchy-2008)
- [Donald B. Rubin (1987) — Multiple Imputation for Nonresponse in Surveys](../CITATIONS.md#rubin-mi-1987)
- [Stef van Buuren and Karin Groothuis-Oudshoorn (2011) — mice: Multivariate Imputation by Chained Equations in R](../CITATIONS.md#vanbuuren-mice-2011)
- [John Barnard and Donald B. Rubin (1999) — Small-sample degrees of freedom with multiple imputation](../CITATIONS.md#barnard-rubin-1999)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
