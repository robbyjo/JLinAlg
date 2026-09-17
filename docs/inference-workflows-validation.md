# Differential, EWAS region, multiple-testing, and imputation validation

Source-build addition, September 2026. The implementation is organized under
`org.jlinalg.differential`, `org.jlinalg.ewas`,
`org.jlinalg.multipletesting`, and `org.jlinalg.imputation`; file-oriented
workflows are in `XwasInferenceCli`.

## Executed numerical gates

The focused Gradle suite contains 14 tests and passed locally on 2026-09-17.
It covers the public APIs, CLI dispatch, output sidecars, reproducible streams,
failure accounting, and the following independent or external oracles:

| Contract | Executed check |
| --- | --- |
| Gaussian empirical Bayes | Coefficients agree with frozen limma to `2e-12`; moderated-t correlation exceeds `0.995` |
| Voom | Log-fold-change correlation with frozen limma/voom exceeds `0.97`; every precision weight is positive |
| Negative binomial | Log2-fold-change correlation exceeds `0.90` against both edgeR and DESeq2; all 14 simulated non-null directions agree |
| NB likelihood mechanics | Majority-invariant size-factor fixture checks the fitted log mean ratio and geometric-mean-one normalization |
| Spatial EWAS | Hand evaluation of all nine exponential-covariance terms reproduces the signed regional Z to `1e-12` |
| BH/IHW | A known BH vector is exact to `1e-15`; changing held-out-fold p-values leaves its learned weights unchanged |
| Hierarchy | Weights sum to one, ancestor gating is enforced, and arbitrary-dependence weighted Bonferroni output is checked |
| MICE | Same seed gives byte-identical completed numeric matrices, observed cells remain unchanged, and all cells are completed |
| Rubin pooling | A three-imputation hand calculation gives within variance 4, between variance 1, and total variance `16/3` |

`src/test/resources/inference-workflows/generate-reference.R` generated the
checked-in continuous/count inputs and package outputs with seed `17092026`.
`versions.tsv` freezes R 4.6.1, Bioconductor 3.23, limma 3.68.5, edgeR 4.10.5,
and DESeq2 1.52.0. No cohort or downloaded analysis data are included.

The local reference run emitted DESeq2's documented fallback from a parametric
to a local dispersion trend for this small synthetic matrix. That behavior is
part of the frozen comparator, not silently repaired in JLinAlg.

## Reproduce

Install reference packages into a disposable repository-local library if they
are not already present:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -e "dir.create('build/r-lib',recursive=TRUE,showWarnings=FALSE); .libPaths(c(normalizePath('build/r-lib'),.libPaths())); if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org'); BiocManager::install(c('limma','edgeR','DESeq2'),ask=FALSE,update=FALSE)"
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/inference-workflows/generate-reference.R
```

Run the focused and full Java gates serially so concurrent Gradle jobs do not
race over `build/`:

```powershell
.\gradlew.bat test --tests "org.jlinalg.differential.*" --tests "org.jlinalg.ewas.*" --tests "org.jlinalg.multipletesting.*" --tests "org.jlinalg.imputation.*" --tests "org.jlinalg.cli.XwasInferenceCliTest" --no-daemon --no-parallel
.\gradlew.bat clean check assemble --no-daemon --no-parallel
```

## Supported scope and non-equivalence boundaries

- `limma` uses a common scaled-inverse-chi-square residual-variance prior.
  Robust/winsorized priors, duplicate correlation, array weights, treat-style
  fold-change testing, and arbitrary contrast matrices in the CLI are outside
  this first interface.
- `voom` uses median-ratio library scaling and a direct local precision trend.
  It does not claim exact TMM or lowess identity with limma/edgeR for every
  workload. The separate NB path must be used for count likelihood inference.
- The NB path uses fixed shrunken dispersions in model-based Wald inference.
  edgeR quasi-likelihood, DESeq2 coefficient priors, Cook outlier replacement,
  independent filtering, and their full dispersion-trend choices are not
  relabeled as implemented.
- EWAS regions are nonoverlapping coordinate runs using a declared exponential
  covariance. The build is required metadata; coordinate lifting and reference
  annotation are upstream responsibilities. This is not comb-p or DMRcate
  option-for-option parity.
- IHW requires a covariate independent of null p-values and independent folds.
  The compact cross-fitted bin estimator is not the complete IHW convex
  optimizer. Hierarchical output controls FWER, not FDR.
- MICE uses predictive mean matching, binary logistic draws, and categorical
  probability draws for numeric or CLI-encoded categories. It does not include
  multilevel random effects, passive formulas, survey weights, or MNAR models.
  Chain summaries are diagnostics, not proof of convergence or correct MAR
  specification.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Gordon K. Smyth (2004) — Linear models and empirical Bayes methods for assessing differential expression in microarray experiments](CITATIONS.md#smyth-limma-2004) — [PMID: 16646809](https://pubmed.ncbi.nlm.nih.gov/16646809/)
- [Charity W. Law, Yunshun Chen, Wei Shi, and Gordon K. Smyth (2014) — voom: precision weights unlock linear model analysis tools for RNA-seq read counts](CITATIONS.md#law-voom-2014) — [PMID: 24485249](https://pubmed.ncbi.nlm.nih.gov/24485249/) · [PMCID: PMC4053721](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053721/)
- [Mark D. Robinson, Davis J. McCarthy, and Gordon K. Smyth (2010) — edgeR: a Bioconductor package for differential expression analysis of digital gene expression data](CITATIONS.md#robinson-edger-2010) — [PMID: 19910308](https://pubmed.ncbi.nlm.nih.gov/19910308/) · [PMCID: PMC2796818](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/)
- [Michael I. Love, Wolfgang Huber, and Simon Anders (2014) — Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](CITATIONS.md#love-deseq2-2014) — [PMID: 25516281](https://pubmed.ncbi.nlm.nih.gov/25516281/) · [PMCID: PMC4302049](https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/)
- [Brent S. Pedersen, David A. Schwartz, Ivana V. Yang, and Katerina J. Kechris (2012) — Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values](CITATIONS.md#pedersen-combp-2012) — [PMID: 22954632](https://pubmed.ncbi.nlm.nih.gov/22954632/) · [PMCID: PMC3496335](https://pmc.ncbi.nlm.nih.gov/articles/PMC3496335/)
- [Timothy J. Peters et al. (2015) — De novo identification of differentially methylated regions in the human genome](CITATIONS.md#peters-dmrcate-2015) — [PMID: 25972926](https://pubmed.ncbi.nlm.nih.gov/25972926/) · [PMCID: PMC4429355](https://pmc.ncbi.nlm.nih.gov/articles/PMC4429355/)
- [Yoav Benjamini and Yosef Hochberg (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](CITATIONS.md#benjamini-hochberg-1995)
- [Nikolaos Ignatiadis, Bernd Klaus, Judith B. Zaugg, and Wolfgang Huber (2016) — Data-driven hypothesis weighting increases detection power in genome-scale multiple testing](CITATIONS.md#ignatiadis-ihw-2016) — [PMID: 27240256](https://pubmed.ncbi.nlm.nih.gov/27240256/) · [PMCID: PMC4930141](https://pmc.ncbi.nlm.nih.gov/articles/PMC4930141/)
- [Nicolai Meinshausen (2008) — Hierarchical testing of variable importance](CITATIONS.md#meinshausen-hierarchy-2008)
- [Donald B. Rubin (1987) — Multiple Imputation for Nonresponse in Surveys](CITATIONS.md#rubin-mi-1987)
- [Stef van Buuren and Karin Groothuis-Oudshoorn (2011) — mice: Multivariate Imputation by Chained Equations in R](CITATIONS.md#vanbuuren-mice-2011)
- [John Barnard and Donald B. Rubin (1999) — Small-sample degrees of freedom with multiple imputation](CITATIONS.md#barnard-rubin-1999)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
