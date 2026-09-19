# Scientific citations

This bibliography links JLinAlg features to the primary publications that introduced or established their statistical methods. Cite both the relevant method paper and JLinAlg when reporting an analysis. Inclusion here documents methodological provenance; it is not a claim that every implementation detail is identical to the cited software.

The website provides a [searchable citation index](https://robbyjo.github.io/JLinAlg/citations.html). This file and the per-vignette citation panels are generated from [citations.json](citations.json).

PMID and PMCID values are retrieved from [NCBI PubMed](https://pubmed.ncbi.nlm.nih.gov/) with [enrich-citation-identifiers.py](enrich-citation-identifiers.py).

## v0.3.6 inference additions

The v0.3.6 release adds four connected inference workflows. These quick links expose their primary methodological foundations; the canonical entries remain in the topic bibliography below.

### Empirical-Bayes differential analysis

Moderated Gaussian models, precision-weighted counts, and negative-binomial dispersion shrinkage.

- [Smyth (2004) — Linear models and empirical Bayes methods for assessing differential expression in microarray experiments](#smyth-limma-2004)
- [Law (2014) — voom: precision weights unlock linear model analysis tools for RNA-seq read counts](#law-voom-2014)
- [Robinson (2010) — edgeR: a Bioconductor package for differential expression analysis of digital gene expression data](#robinson-edger-2010)
- [Love (2014) — Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](#love-deseq2-2014)

### Region-level EWAS

Coordinate-aware aggregation with explicit spatial dependence and region-wide multiplicity control.

- [Pedersen (2012) — Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values](#pedersen-combp-2012)
- [Peters et al. (2015) — De novo identification of differentially methylated regions in the human genome](#peters-dmrcate-2015)

### Adaptive and hierarchical testing

Weighted FDR and prespecified hierarchy-aware familywise inference.

- [Benjamini (1995) — Controlling the false discovery rate: a practical and powerful approach to multiple testing](#benjamini-hochberg-1995)
- [Ignatiadis (2016) — Data-driven hypothesis weighting increases detection power in genome-scale multiple testing](#ignatiadis-ihw-2016)
- [Meinshausen (2008) — Hierarchical testing of variable importance](#meinshausen-hierarchy-2008)

### Multiple-imputation inference

Chained equations, Rubin variance pooling, and finite-sample degrees of freedom.

- [Rubin (1987) — Multiple Imputation for Nonresponse in Surveys](#rubin-mi-1987)
- [Buuren (2011) — mice: Multivariate Imputation by Chained Equations in R](#vanbuuren-mice-2011)
- [Barnard (1999) — Small-sample degrees of freedom with multiple imputation](#barnard-rubin-1999)

## Foundational inference

<a id="student-1908"></a>
- **Student [W. S. Gosset] (1908).** The probable error of a mean. *Biometrika 6:1-25.* [DOI: 10.2307/2331554](https://doi.org/10.2307/2331554)
  - JLinAlg methods: Student t test, small-sample mean inference.

<a id="welch-1947"></a>
- **B. L. Welch (1947).** The generalization of Student's problem when several different population variances are involved. *Biometrika 34:28-35.* [DOI: 10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28) · [PMID: 20287819](https://pubmed.ncbi.nlm.nih.gov/20287819/)
  - JLinAlg methods: Welch t test, unequal-variance inference.

<a id="pearson-1895"></a>
- **Karl Pearson (1895).** Notes on regression and inheritance in the case of two parents. *Proceedings of the Royal Society of London 58:240-242.* [DOI: 10.1098/rspl.1895.0041](https://doi.org/10.1098/rspl.1895.0041)
  - JLinAlg methods: Pearson correlation, linear association.

<a id="spearman-1904"></a>
- **Charles Spearman (1904).** The proof and measurement of association between two things. *American Journal of Psychology 15:72-101.* [DOI: 10.2307/1412159](https://doi.org/10.2307/1412159)
  - JLinAlg methods: Spearman correlation, rank association.

<a id="wilcoxon-1945"></a>
- **Frank Wilcoxon (1945).** Individual comparisons by ranking methods. *Biometrics Bulletin 1:80-83.* [DOI: 10.2307/3001968](https://doi.org/10.2307/3001968)
  - JLinAlg methods: Wilcoxon signed-rank test, rank inference.

<a id="mann-whitney-1947"></a>
- **H. B. Mann and D. R. Whitney (1947).** On a test of whether one of two random variables is stochastically larger than the other. *Annals of Mathematical Statistics 18:50-60.* [DOI: 10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)
  - JLinAlg methods: Mann-Whitney U test, two-sample rank inference.

<a id="fisher-1922"></a>
- **R. A. Fisher (1922).** On the interpretation of chi-square from contingency tables, and the calculation of P. *Journal of the Royal Statistical Society 85:87-94.* [DOI: 10.2307/2340521](https://doi.org/10.2307/2340521)
  - JLinAlg methods: Fisher exact test, contingency tables.

<a id="mcnemar-1947"></a>
- **Quinn McNemar (1947).** Note on the sampling error of the difference between correlated proportions or percentages. *Psychometrika 12:153-157.* [DOI: 10.1007/BF02295996](https://doi.org/10.1007/BF02295996) · [PMID: 20254758](https://pubmed.ncbi.nlm.nih.gov/20254758/)
  - JLinAlg methods: McNemar test, paired binary data.

<a id="friedman-1937"></a>
- **Milton Friedman (1937).** The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *Journal of the American Statistical Association 32:675-701.* [DOI: 10.1080/01621459.1937.10503522](https://doi.org/10.1080/01621459.1937.10503522)
  - JLinAlg methods: Friedman test, blocked rank inference.

<a id="benjamini-hochberg-1995"></a>
- **Yoav Benjamini and Yosef Hochberg (1995).** Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B 57:289-300.* [DOI: 10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
  - JLinAlg methods: Benjamini-Hochberg FDR, multiple testing.

## Regression

<a id="nelder-wedderburn-1972"></a>
- **J. A. Nelder and R. W. M. Wedderburn (1972).** Generalized linear models. *Journal of the Royal Statistical Society A 135:370-384.* [DOI: 10.2307/2344614](https://doi.org/10.2307/2344614)
  - JLinAlg methods: generalized linear models, IRLS, logistic regression, Poisson regression.

<a id="bliss-1934"></a>
- **C. I. Bliss (1934).** The method of probits. *Science 79:38-39.* [DOI: 10.1126/science.79.2037.38](https://doi.org/10.1126/science.79.2037.38) · [PMID: 17813446](https://pubmed.ncbi.nlm.nih.gov/17813446/)
  - JLinAlg methods: probit regression, binary response.

<a id="hoerl-kennard-1970"></a>
- **Arthur E. Hoerl and Robert W. Kennard (1970).** Ridge regression: biased estimation for nonorthogonal problems. *Technometrics 12:55-67.* [DOI: 10.1080/00401706.1970.10488634](https://doi.org/10.1080/00401706.1970.10488634)
  - JLinAlg methods: ridge regression, prediction scores.

<a id="tibshirani-1996"></a>
- **Robert Tibshirani (1996).** Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society B 58:267-288.* [DOI: 10.1111/j.2517-6161.1996.tb02080.x](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x)
  - JLinAlg methods: lasso, sparse regression.

<a id="zou-hastie-2005"></a>
- **Hui Zou and Trevor Hastie (2005).** Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society B 67:301-320.* [DOI: 10.1111/j.1467-9868.2005.00503.x](https://doi.org/10.1111/j.1467-9868.2005.00503.x)
  - JLinAlg methods: elastic net, grouped shrinkage.

<a id="lee-2016"></a>
- **Jason D. Lee, Dennis L. Sun, Yuekai Sun, and Jonathan E. Taylor (2016).** Exact post-selection inference, with application to the lasso. *Annals of Statistics 44:907-927.* [DOI: 10.1214/15-AOS1371](https://doi.org/10.1214/15-AOS1371)
  - JLinAlg methods: polyhedral selective inference, lasso inference.

<a id="koenker-bassett-1978"></a>
- **Roger Koenker and Gilbert Bassett Jr. (1978).** Regression quantiles. *Econometrica 46:33-50.* [DOI: 10.2307/1913643](https://doi.org/10.2307/1913643)
  - JLinAlg methods: quantile regression.

<a id="friedman-1984"></a>
- **Jerome H. Friedman (1984).** A variable span smoother. *Laboratory for Computational Statistics Technical Report 5.* [DOI: 10.2172/1447470](https://doi.org/10.2172/1447470)
  - JLinAlg methods: supersmoother, adaptive span smoothing.

<a id="nadaraya-1964"></a>
- **E. A. Nadaraya (1964).** On estimating regression. *Theory of Probability and Its Applications 9:141-142.* [DOI: 10.1137/1109020](https://doi.org/10.1137/1109020)
  - JLinAlg methods: kernel regression, Nadaraya-Watson estimator.

<a id="robinson-1988"></a>
- **Peter M. Robinson (1988).** Root-N-consistent semiparametric regression. *Econometrica 56:931-954.* [DOI: 10.2307/1912705](https://doi.org/10.2307/1912705)
  - JLinAlg methods: partially linear regression, semiparametric regression.

<a id="cleveland-1979"></a>
- **William S. Cleveland (1979).** Robust locally weighted regression and smoothing scatterplots. *Journal of the American Statistical Association 74:829-836.* [DOI: 10.1080/01621459.1979.10481038](https://doi.org/10.1080/01621459.1979.10481038)
  - JLinAlg methods: LOESS, robust local regression.

<a id="hastie-tibshirani-1986"></a>
- **Trevor Hastie and Robert Tibshirani (1986).** Generalized additive models. *Statistical Science 1:297-310.* [DOI: 10.1214/ss/1177013604](https://doi.org/10.1214/ss/1177013604)
  - JLinAlg methods: generalized additive models, additive smooths.

<a id="eilers-marx-1996"></a>
- **Paul H. C. Eilers and Brian D. Marx (1996).** Flexible smoothing with B-splines and penalties. *Statistical Science 11:89-121.* [DOI: 10.1214/ss/1038425655](https://doi.org/10.1214/ss/1038425655)
  - JLinAlg methods: P-splines, penalized B-splines.

<a id="wood-2011"></a>
- **Simon N. Wood (2011).** Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *Journal of the Royal Statistical Society B 73:3-36.* [DOI: 10.1111/j.1467-9868.2010.00749.x](https://doi.org/10.1111/j.1467-9868.2010.00749.x)
  - JLinAlg methods: GAM REML, smoothness selection, penalized GLM.

<a id="rigby-stasinopoulos-2005"></a>
- **Robert A. Rigby and D. Mikis Stasinopoulos (2005).** Generalized additive models for location, scale and shape. *Applied Statistics 54:507-554.* [DOI: 10.1111/j.1467-9876.2005.00510.x](https://doi.org/10.1111/j.1467-9876.2005.00510.x)
  - JLinAlg methods: GAMLSS, distributional regression.

<a id="yee-2010"></a>
- **Thomas W. Yee (2010).** The VGAM package for categorical data analysis. *Journal of Statistical Software 32:1-34.* [DOI: 10.18637/jss.v032.i10](https://doi.org/10.18637/jss.v032.i10)
  - JLinAlg methods: VGAM, ordinal and categorical regression.

<a id="ferrari-cribari-neto-2004"></a>
- **Silvia Ferrari and Francisco Cribari-Neto (2004).** Beta regression for modelling rates and proportions. *Journal of Applied Statistics 31:799-815.* [DOI: 10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)
  - JLinAlg methods: beta regression, mean-precision models.

<a id="beasley-2009"></a>
- **T. Mark Beasley, Stephen Erickson, and David B. Allison (2009).** Rank-based inverse normal transformations are increasingly used, but are they merited?. *Behavior Genetics 39:580-595.* [DOI: 10.1007/s10519-009-9281-0](https://doi.org/10.1007/s10519-009-9281-0) · [PMID: 19526352](https://pubmed.ncbi.nlm.nih.gov/19526352/) · [PMCID: PMC2921808](https://pmc.ncbi.nlm.nih.gov/articles/PMC2921808/)
  - JLinAlg methods: rank inverse-normal transform, omics transforms.

## Mixed, longitudinal, and survival models

<a id="patterson-thompson-1971"></a>
- **H. D. Patterson and Robin Thompson (1971).** Recovery of inter-block information when block sizes are unequal. *Biometrika 58:545-554.* [DOI: 10.1093/biomet/58.3.545](https://doi.org/10.1093/biomet/58.3.545)
  - JLinAlg methods: REML, variance components.

<a id="henderson-1975"></a>
- **C. R. Henderson (1975).** Best linear unbiased estimation and prediction under a selection model. *Biometrics 31:423-447.* [DOI: 10.2307/2529430](https://doi.org/10.2307/2529430)
  - JLinAlg methods: BLUP, mixed-model equations.

<a id="kenward-roger-1997"></a>
- **Michael G. Kenward and James H. Roger (1997).** Small sample inference for fixed effects from restricted maximum likelihood. *Biometrics 53:983-997.* [DOI: 10.2307/2533558](https://doi.org/10.2307/2533558)
  - JLinAlg methods: Kenward-Roger inference, small-sample mixed models.

<a id="breslow-clayton-1993"></a>
- **Norman E. Breslow and David G. Clayton (1993).** Approximate inference in generalized linear mixed models. *Journal of the American Statistical Association 88:9-25.* [DOI: 10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284)
  - JLinAlg methods: PQL, generalized linear mixed models.

<a id="lambert-1992"></a>
- **Diane Lambert (1992).** Zero-inflated Poisson regression, with an application to defects in manufacturing. *Technometrics 34:1-14.* [DOI: 10.1080/00401706.1992.10485228](https://doi.org/10.1080/00401706.1992.10485228)
  - JLinAlg methods: zero-inflated Poisson, mixture count models.

<a id="lindstrom-bates-1990"></a>
- **Mary J. Lindstrom and Douglas M. Bates (1990).** Nonlinear mixed effects models for repeated measures data. *Biometrics 46:673-687.* [DOI: 10.2307/2532087](https://doi.org/10.2307/2532087)
  - JLinAlg methods: nonlinear mixed-effects models.

<a id="liang-zeger-1986"></a>
- **Kung-Yee Liang and Scott L. Zeger (1986).** Longitudinal data analysis using generalized linear models. *Biometrika 73:13-22.* [DOI: 10.1093/biomet/73.1.13](https://doi.org/10.1093/biomet/73.1.13)
  - JLinAlg methods: GEE, working correlation, sandwich inference.

<a id="cox-1972"></a>
- **D. R. Cox (1972).** Regression models and life-tables. *Journal of the Royal Statistical Society B 34:187-220.* [DOI: 10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)
  - JLinAlg methods: Cox proportional hazards, partial likelihood.

<a id="efron-1977"></a>
- **Bradley Efron (1977).** The efficiency of Cox's likelihood function for censored data. *Journal of the American Statistical Association 72:557-565.* [DOI: 10.1080/01621459.1977.10480613](https://doi.org/10.1080/01621459.1977.10480613)
  - JLinAlg methods: Efron ties, Cox regression.

## Quantitative and statistical genetics

<a id="henderson-1976"></a>
- **C. R. Henderson (1976).** A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. *Biometrics 32:69-83.* [DOI: 10.2307/2529339](https://doi.org/10.2307/2529339)
  - JLinAlg methods: pedigree A inverse, animal models.

<a id="vanraden-2008"></a>
- **Paul M. VanRaden (2008).** Efficient methods to compute genomic predictions. *Journal of Dairy Science 91:4414-4423.* [DOI: 10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980) · [PMID: 18946147](https://pubmed.ncbi.nlm.nih.gov/18946147/)
  - JLinAlg methods: genomic relationship matrix, genomic prediction.

<a id="yang-gcta-2011"></a>
- **Jian Yang et al. (2011).** GCTA: a tool for genome-wide complex trait analysis. *American Journal of Human Genetics 88:76-82.* [DOI: 10.1016/j.ajhg.2010.11.011](https://doi.org/10.1016/j.ajhg.2010.11.011) · [PMID: 21167468](https://pubmed.ncbi.nlm.nih.gov/21167468/) · [PMCID: PMC3014363](https://pmc.ncbi.nlm.nih.gov/articles/PMC3014363/)
  - JLinAlg methods: GRM, SNP heritability, GCTA.

<a id="kang-emmax-2010"></a>
- **Hyun Min Kang et al. (2010).** Variance component model to account for sample structure in genome-wide association studies. *Nature Genetics 42:348-354.* [DOI: 10.1038/ng.548](https://doi.org/10.1038/ng.548) · [PMID: 20208533](https://pubmed.ncbi.nlm.nih.gov/20208533/) · [PMCID: PMC3092069](https://pmc.ncbi.nlm.nih.gov/articles/PMC3092069/)
  - JLinAlg methods: EMMAX, mixed-model GWAS.

<a id="zhang-p3d-2010"></a>
- **Zhiwu Zhang et al. (2010).** Mixed linear model approach adapted for genome-wide association studies. *Nature Genetics 42:355-360.* [DOI: 10.1038/ng.546](https://doi.org/10.1038/ng.546) · [PMID: 20208535](https://pubmed.ncbi.nlm.nih.gov/20208535/) · [PMCID: PMC2931336](https://pmc.ncbi.nlm.nih.gov/articles/PMC2931336/)
  - JLinAlg methods: P3D, compressed mixed-model GWAS.

<a id="wu-skat-2011"></a>
- **Michael C. Wu et al. (2011).** Rare-variant association testing for sequencing data with the sequence kernel association test. *American Journal of Human Genetics 89:82-93.* [DOI: 10.1016/j.ajhg.2011.05.029](https://doi.org/10.1016/j.ajhg.2011.05.029) · [PMID: 21737059](https://pubmed.ncbi.nlm.nih.gov/21737059/) · [PMCID: PMC3135811](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135811/)
  - JLinAlg methods: SKAT, rare-variant set tests.

<a id="lee-skato-2012"></a>
- **Seunggeun Lee et al. (2012).** Optimal unified approach for rare-variant association testing. *American Journal of Human Genetics 91:224-237.* [DOI: 10.1016/j.ajhg.2012.06.007](https://doi.org/10.1016/j.ajhg.2012.06.007) · [PMID: 22863193](https://pubmed.ncbi.nlm.nih.gov/22863193/) · [PMCID: PMC3415556](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
  - JLinAlg methods: SKAT-O, burden-kernel omnibus tests.

<a id="liu-acat-2019"></a>
- **Yaowu Liu et al. (2019).** ACAT: A fast and powerful p value combination method for rare-variant analysis in sequencing studies. *American Journal of Human Genetics 104:410-421.* [DOI: 10.1016/j.ajhg.2019.01.002](https://doi.org/10.1016/j.ajhg.2019.01.002) · [PMID: 30849328](https://pubmed.ncbi.nlm.nih.gov/30849328/) · [PMCID: PMC6407498](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/)
  - JLinAlg methods: ACAT, ACAT-V, ACAT-O, rare-variant omnibus tests.

<a id="madsen-browning-2009"></a>
- **B. E. Madsen and S. R. Browning (2009).** A groupwise association test for rare mutations using a weighted sum statistic. *PLoS Genetics 5:e1000384.* [DOI: 10.1371/journal.pgen.1000384](https://doi.org/10.1371/journal.pgen.1000384) · [PMID: 19214210](https://pubmed.ncbi.nlm.nih.gov/19214210/) · [PMCID: PMC2633048](https://pmc.ncbi.nlm.nih.gov/articles/PMC2633048/)
  - JLinAlg methods: weighted burden test, rare variants.

<a id="yang-cojo-2012"></a>
- **Jian Yang et al. (2012).** Conditional and joint multiple-SNP analysis of GWAS summary statistics. *Nature Genetics 44:369-375.* [DOI: 10.1038/ng.2213](https://doi.org/10.1038/ng.2213) · [PMID: 22426310](https://pubmed.ncbi.nlm.nih.gov/22426310/) · [PMCID: PMC3593158](https://pmc.ncbi.nlm.nih.gov/articles/PMC3593158/)
  - JLinAlg methods: COJO, conditional summary statistics.

<a id="bulik-sullivan-ldsc-2015"></a>
- **Brendan K. Bulik-Sullivan et al. (2015).** LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics 47:291-295.* [DOI: 10.1038/ng.3211](https://doi.org/10.1038/ng.3211) · [PMID: 25642630](https://pubmed.ncbi.nlm.nih.gov/25642630/) · [PMCID: PMC4495769](https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/)
  - JLinAlg methods: LD Score regression, SNP heritability, confounding intercept.

<a id="bulik-sullivan-rg-2015"></a>
- **Brendan K. Bulik-Sullivan et al. (2015).** An atlas of genetic correlations across human diseases and traits. *Nature Genetics 47:1236-1241.* [DOI: 10.1038/ng.3406](https://doi.org/10.1038/ng.3406) · [PMID: 26414676](https://pubmed.ncbi.nlm.nih.gov/26414676/) · [PMCID: PMC4797329](https://pmc.ncbi.nlm.nih.gov/articles/PMC4797329/)
  - JLinAlg methods: genetic correlation, bivariate LDSC.

<a id="gamazon-predixcan-2015"></a>
- **Eric R. Gamazon et al. (2015).** A gene-based association method for mapping traits using reference transcriptome data. *Nature Genetics 47:1091-1098.* [DOI: 10.1038/ng.3367](https://doi.org/10.1038/ng.3367) · [PMID: 26258848](https://pubmed.ncbi.nlm.nih.gov/26258848/) · [PMCID: PMC4552594](https://pmc.ncbi.nlm.nih.gov/articles/PMC4552594/)
  - JLinAlg methods: PrediXcan, genetically predicted expression, TWAS.

<a id="gusev-fusion-2016"></a>
- **Alexander Gusev et al. (2016).** Integrative approaches for large-scale transcriptome-wide association studies. *Nature Genetics 48:245-252.* [DOI: 10.1038/ng.3506](https://doi.org/10.1038/ng.3506) · [PMID: 26854917](https://pubmed.ncbi.nlm.nih.gov/26854917/) · [PMCID: PMC4767558](https://pmc.ncbi.nlm.nih.gov/articles/PMC4767558/)
  - JLinAlg methods: FUSION TWAS, summary TWAS, PWAS.

<a id="grotzinger-genomic-sem-2019"></a>
- **Andrew D. Grotzinger et al. (2019).** Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nature Human Behaviour 3:513-525.* [DOI: 10.1038/s41562-019-0566-x](https://doi.org/10.1038/s41562-019-0566-x) · [PMID: 30962613](https://pubmed.ncbi.nlm.nih.gov/30962613/) · [PMCID: PMC6520146](https://pmc.ncbi.nlm.nih.gov/articles/PMC6520146/)
  - JLinAlg methods: Genomic SEM, shared genetic factors.

<a id="euesden-prsice-2015"></a>
- **Jack Euesden, Cathryn M. Lewis, and Paul F. O'Reilly (2015).** PRSice: polygenic risk score software. *Bioinformatics 31:1466-1468.* [DOI: 10.1093/bioinformatics/btu848](https://doi.org/10.1093/bioinformatics/btu848) · [PMID: 25550326](https://pubmed.ncbi.nlm.nih.gov/25550326/) · [PMCID: PMC4410663](https://pmc.ncbi.nlm.nih.gov/articles/PMC4410663/)
  - JLinAlg methods: polygenic scores, score evaluation.

## Latent variation and batch correction

<a id="leek-storey-2007"></a>
- **Jeffrey T. Leek and John D. Storey (2007).** Capturing heterogeneity in gene expression studies by surrogate variable analysis. *PLoS Genetics 3:e161.* [DOI: 10.1371/journal.pgen.0030161](https://doi.org/10.1371/journal.pgen.0030161) · [PMID: 17907809](https://pubmed.ncbi.nlm.nih.gov/17907809/) · [PMCID: PMC1994707](https://pmc.ncbi.nlm.nih.gov/articles/PMC1994707/)
  - JLinAlg methods: SVA, latent confounders.

<a id="stegle-peer-2012"></a>
- **Oliver Stegle et al. (2012).** Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. *Nature Protocols 7:500-507.* [DOI: 10.1038/nprot.2011.457](https://doi.org/10.1038/nprot.2011.457) · [PMID: 22343431](https://pubmed.ncbi.nlm.nih.gov/22343431/) · [PMCID: PMC3398141](https://pmc.ncbi.nlm.nih.gov/articles/PMC3398141/)
  - JLinAlg methods: PEER, Bayesian factor analysis.

<a id="johnson-combat-2007"></a>
- **W. Evan Johnson, Cheng Li, and Ariel Rabinovic (2007).** Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics 8:118-127.* [DOI: 10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037) · [PMID: 16632515](https://pubmed.ncbi.nlm.nih.gov/16632515/)
  - JLinAlg methods: ComBat, batch correction, empirical Bayes.

## Causal and mediation analysis

<a id="baron-kenny-1986"></a>
- **Reuben M. Baron and David A. Kenny (1986).** The moderator-mediator variable distinction in social psychological research. *Journal of Personality and Social Psychology 51:1173-1182.* [DOI: 10.1037/0022-3514.51.6.1173](https://doi.org/10.1037/0022-3514.51.6.1173)
  - JLinAlg methods: mediation, path decomposition.

<a id="sobel-1982"></a>
- **Michael E. Sobel (1982).** Asymptotic confidence intervals for indirect effects in structural equation models. *Sociological Methodology 13:290-312.* [DOI: 10.2307/270723](https://doi.org/10.2307/270723)
  - JLinAlg methods: Sobel test, indirect-effect inference.

<a id="staiger-stock-1997"></a>
- **Douglas Staiger and James H. Stock (1997).** Instrumental variables regression with weak instruments. *Econometrica 65:557-586.* [DOI: 10.2307/2171753](https://doi.org/10.2307/2171753)
  - JLinAlg methods: 2SLS, weak instruments, first-stage diagnostics.

<a id="white-1980"></a>
- **Halbert White (1980).** A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica 48:817-838.* [DOI: 10.2307/1912934](https://doi.org/10.2307/1912934)
  - JLinAlg methods: heteroskedasticity-robust covariance, robust 2SLS inference.

<a id="bowden-mr-egger-2015"></a>
- **Jack Bowden, George Davey Smith, and Stephen Burgess (2015).** Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *International Journal of Epidemiology 44:512-525.* [DOI: 10.1093/ije/dyv080](https://doi.org/10.1093/ije/dyv080) · [PMID: 26050253](https://pubmed.ncbi.nlm.nih.gov/26050253/) · [PMCID: PMC4469799](https://pmc.ncbi.nlm.nih.gov/articles/PMC4469799/)
  - JLinAlg methods: MR-Egger, directional pleiotropy, Egger intercept.

<a id="bowden-weighted-median-2016"></a>
- **Jack Bowden et al. (2016).** Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator. *Genetic Epidemiology 40:304-314.* [DOI: 10.1002/gepi.21965](https://doi.org/10.1002/gepi.21965) · [PMID: 27061298](https://pubmed.ncbi.nlm.nih.gov/27061298/) · [PMCID: PMC4849733](https://pmc.ncbi.nlm.nih.gov/articles/PMC4849733/)
  - JLinAlg methods: weighted-median MR, invalid instruments.

<a id="hemani-steiger-2017"></a>
- **Gibran Hemani, Kate Tilling, and George Davey Smith (2017).** Orienting the causal relationship between imprecisely measured traits using GWAS summary data. *PLoS Genetics 13:e1007081.* [DOI: 10.1371/journal.pgen.1007081](https://doi.org/10.1371/journal.pgen.1007081) · [PMID: 29149188](https://pubmed.ncbi.nlm.nih.gov/29149188/) · [PMCID: PMC5711033](https://pmc.ncbi.nlm.nih.gov/articles/PMC5711033/)
  - JLinAlg methods: Steiger directionality, MR orientation.

<a id="zhao-mr-raps-2020"></a>
- **Qingyuan Zhao et al. (2020).** Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score. *Annals of Statistics 48:1742-1769.* [DOI: 10.1214/19-AOS1866](https://doi.org/10.1214/19-AOS1866)
  - JLinAlg methods: MR-RAPS, weak instruments, idiosyncratic pleiotropy.

<a id="verbanck-mr-presso-2018"></a>
- **Marie Verbanck et al. (2018).** Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization. *Nature Genetics 50:693-698.* [DOI: 10.1038/s41588-018-0099-7](https://doi.org/10.1038/s41588-018-0099-7) · [PMID: 29686387](https://pubmed.ncbi.nlm.nih.gov/29686387/) · [PMCID: PMC6083837](https://pmc.ncbi.nlm.nih.gov/articles/PMC6083837/)
  - JLinAlg methods: MR-PRESSO, pleiotropic outliers.

<a id="burgess-mvmr-2015"></a>
- **Stephen Burgess and Simon G. Thompson (2015).** Multivariable Mendelian randomization: the use of pleiotropic genetic variants to estimate causal effects. *American Journal of Epidemiology 181:251-260.* [DOI: 10.1093/aje/kwu283](https://doi.org/10.1093/aje/kwu283) · [PMID: 25632051](https://pubmed.ncbi.nlm.nih.gov/25632051/) · [PMCID: PMC4325677](https://pmc.ncbi.nlm.nih.gov/articles/PMC4325677/)
  - JLinAlg methods: multivariable MR, direct causal effects.

<a id="sanderson-mvmr-diagnostics-2021"></a>
- **Eleanor Sanderson, Wesley Spiller, and Jack Bowden (2021).** Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization. *Statistics in Medicine 40:5434-5452.* [DOI: 10.1002/sim.9133](https://doi.org/10.1002/sim.9133) · [PMID: 34338327](https://pubmed.ncbi.nlm.nih.gov/34338327/) · [PMCID: PMC9479726](https://pmc.ncbi.nlm.nih.gov/articles/PMC9479726/)
  - JLinAlg methods: multivariable MR diagnostics, conditional instrument strength, pleiotropy heterogeneity.

<a id="zuber-mr2-2023"></a>
- **Verena Zuber et al. (2023).** Multi-response Mendelian randomization: Identification of shared and distinct exposures for multimorbidity and multiple related disease outcomes. *American Journal of Human Genetics 110:1177-1199.* [DOI: 10.1016/j.ajhg.2023.06.005](https://doi.org/10.1016/j.ajhg.2023.06.005) · [PMID: 37419091](https://pubmed.ncbi.nlm.nih.gov/37419091/) · [PMCID: PMC10357504](https://pmc.ncbi.nlm.nih.gov/articles/PMC10357504/)
  - JLinAlg methods: multi-response MR literature, sparse Bayesian Gaussian-copula MR2 method boundary.

<a id="deng-mrmo-2023"></a>
- **Yangqing Deng et al. (2023).** Two-stage multivariate Mendelian randomization on multiple outcomes with mixed distributions. *Statistical Methods in Medical Research 32:1543-1558.* [DOI: 10.1177/09622802231181220](https://doi.org/10.1177/09622802231181220) · [PMID: 37338962](https://pubmed.ncbi.nlm.nih.gov/37338962/) · [PMCID: PMC10515454](https://pmc.ncbi.nlm.nih.gov/articles/PMC10515454/)
  - JLinAlg methods: individual-level mixed-response MR literature, MRMO method boundary.

<a id="zhang-multivariate-mr-2026"></a>
- **Yuankai Zhang et al. (2026).** Multivariate Mendelian randomization for joint inferences of correlated outcomes. *European Journal of Epidemiology.* [DOI: 10.1007/s10654-026-01406-1](https://doi.org/10.1007/s10654-026-01406-1) · [PMID: 42207415](https://pubmed.ncbi.nlm.nih.gov/42207415/)
  - JLinAlg methods: multivariate MR-IVW, joint correlated-outcome inference, multivariate MR-PRESSO.

<a id="burgess-overlap-2016"></a>
- **Stephen Burgess, Neil M. Davies, and Simon G. Thompson (2016).** Bias due to participant overlap in two-sample Mendelian randomization. *Genetic Epidemiology 40:597-608.* [DOI: 10.1002/gepi.21998](https://doi.org/10.1002/gepi.21998) · [PMID: 27625185](https://pubmed.ncbi.nlm.nih.gov/27625185/) · [PMCID: PMC5082560](https://pmc.ncbi.nlm.nih.gov/articles/PMC5082560/)
  - JLinAlg methods: sample overlap, two-sample MR.

<a id="hemani-mrbase-2018"></a>
- **Gibran Hemani et al. (2018).** The MR-Base platform supports systematic causal inference across the human phenome. *eLife 7:e34408.* [DOI: 10.7554/eLife.34408](https://doi.org/10.7554/eLife.34408) · [PMID: 29846171](https://pubmed.ncbi.nlm.nih.gov/29846171/) · [PMCID: PMC5976434](https://pmc.ncbi.nlm.nih.gov/articles/PMC5976434/)
  - JLinAlg methods: TwoSampleMR, MR-Base, phenome-wide MR.

## Evidence synthesis

<a id="dersimonian-laird-1986"></a>
- **Rebecca DerSimonian and Nan Laird (1986).** Meta-analysis in clinical trials. *Controlled Clinical Trials 7:177-188.* [DOI: 10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456(86)90046-2) · [PMID: 3802833](https://pubmed.ncbi.nlm.nih.gov/3802833/)
  - JLinAlg methods: random-effects meta-analysis, DerSimonian-Laird.

<a id="paule-mandel-1982"></a>
- **Robert C. Paule and John Mandel (1982).** Consensus values and weighting factors. *Journal of Research of the National Bureau of Standards 87:377-385.* [DOI: 10.6028/jres.087.022](https://doi.org/10.6028/jres.087.022) · [PMID: 34566088](https://pubmed.ncbi.nlm.nih.gov/34566088/) · [PMCID: PMC6768160](https://pmc.ncbi.nlm.nih.gov/articles/PMC6768160/)
  - JLinAlg methods: Paule-Mandel tau squared, random-effects meta-analysis.

<a id="knapp-hartung-2003"></a>
- **Guido Knapp and Joachim Hartung (2003).** Improved tests for a random effects meta-regression with a single covariate. *Statistics in Medicine 22:2693-2710.* [DOI: 10.1002/sim.1482](https://doi.org/10.1002/sim.1482) · [PMID: 12939780](https://pubmed.ncbi.nlm.nih.gov/12939780/)
  - JLinAlg methods: Knapp-Hartung inference, meta-regression.

<a id="higgins-thompson-2002"></a>
- **Julian P. T. Higgins and Simon G. Thompson (2002).** Quantifying heterogeneity in a meta-analysis. *Statistics in Medicine 21:1539-1558.* [DOI: 10.1002/sim.1186](https://doi.org/10.1002/sim.1186) · [PMID: 12111919](https://pubmed.ncbi.nlm.nih.gov/12111919/)
  - JLinAlg methods: I-squared, heterogeneity.

## Time series

<a id="kalman-1960"></a>
- **Rudolf E. Kalman (1960).** A new approach to linear filtering and prediction problems. *Journal of Basic Engineering 82:35-45.* [DOI: 10.1115/1.3662552](https://doi.org/10.1115/1.3662552)
  - JLinAlg methods: Kalman filter, state-space models, diffuse smoothing.

<a id="ljung-box-1978"></a>
- **Greta M. Ljung and George E. P. Box (1978).** On a measure of lack of fit in time series models. *Biometrika 65:297-303.* [DOI: 10.1093/biomet/65.2.297](https://doi.org/10.1093/biomet/65.2.297)
  - JLinAlg methods: Ljung-Box test, residual diagnostics.

<a id="hyndman-khandakar-2008"></a>
- **Rob J. Hyndman and Yeasmin Khandakar (2008).** Automatic time series forecasting: the forecast package for R. *Journal of Statistical Software 27:1-22.* [DOI: 10.18637/jss.v027.i03](https://doi.org/10.18637/jss.v027.i03)
  - JLinAlg methods: automatic ARIMA, stepwise order selection.

## Fine mapping, colocalization, and SEM

<a id="wang-susie-2020"></a>
- **Gao Wang, Abhishek Sarkar, Peter Carbonetto, and Matthew Stephens (2020).** A simple new approach to variable selection in regression, with application to genetic fine mapping. *Journal of the Royal Statistical Society B 82:1273-1300.* [DOI: 10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388) · [PMID: 37220626](https://pubmed.ncbi.nlm.nih.gov/37220626/) · [PMCID: PMC10201948](https://pmc.ncbi.nlm.nih.gov/articles/PMC10201948/)
  - JLinAlg methods: SuSiE, credible sets, fine mapping.

<a id="giambartolomei-coloc-2014"></a>
- **Claudia Giambartolomei et al. (2014).** Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genetics 10:e1004383.* [DOI: 10.1371/journal.pgen.1004383](https://doi.org/10.1371/journal.pgen.1004383) · [PMID: 24830394](https://pubmed.ncbi.nlm.nih.gov/24830394/) · [PMCID: PMC4022491](https://pmc.ncbi.nlm.nih.gov/articles/PMC4022491/)
  - JLinAlg methods: coloc, H0-H4 hypotheses.

<a id="wallace-coloc-2021"></a>
- **Chris Wallace (2021).** A more accurate method for colocalisation analysis allowing for multiple causal variants. *PLoS Genetics 17:e1009440.* [DOI: 10.1371/journal.pgen.1009440](https://doi.org/10.1371/journal.pgen.1009440) · [PMID: 34587156](https://pubmed.ncbi.nlm.nih.gov/34587156/) · [PMCID: PMC8504726](https://pmc.ncbi.nlm.nih.gov/articles/PMC8504726/)
  - JLinAlg methods: SuSiE colocalization, multiple causal variants.

<a id="wright-1921"></a>
- **Sewall Wright (1921).** Correlation and causation. *Journal of Agricultural Research 20:557-585.* [Primary source](https://naldc.nal.usda.gov/catalog/IND43966364)
  - JLinAlg methods: path analysis, structural equation models.

<a id="mcardle-mcdonald-1984"></a>
- **John J. McArdle and Roderick P. McDonald (1984).** Some algebraic properties of the Reticular Action Model for moment structures. *British Journal of Mathematical and Statistical Psychology 37:234-251.* [DOI: 10.1111/j.2044-8317.1984.tb00802.x](https://doi.org/10.1111/j.2044-8317.1984.tb00802.x) · [PMID: 6509005](https://pubmed.ncbi.nlm.nih.gov/6509005/)
  - JLinAlg methods: RAM formulation, SEM covariance structure.

<a id="enders-bandalos-2001"></a>
- **Craig K. Enders and Deborah L. Bandalos (2001).** The relative performance of full information maximum likelihood estimation for missing data in structural equation models. *Structural Equation Modeling 8:430-457.* [DOI: 10.1207/S15328007SEM0803_5](https://doi.org/10.1207/S15328007SEM0803_5)
  - JLinAlg methods: FIML, SEM missing data.

<a id="katsikatsou-2012"></a>
- **Myrsini Katsikatsou et al. (2012).** Pairwise likelihood estimation for factor analysis models with ordinal data. *Computational Statistics and Data Analysis 56:4243-4258.* [DOI: 10.1016/j.csda.2012.04.010](https://doi.org/10.1016/j.csda.2012.04.010)
  - JLinAlg methods: ordinal pairwise likelihood, ordinal SEM.

## Enrichment and annotation

<a id="subramanian-gsea-2005"></a>
- **Aravind Subramanian et al. (2005).** Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences 102:15545-15550.* [DOI: 10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102) · [PMID: 16199517](https://pubmed.ncbi.nlm.nih.gov/16199517/) · [PMCID: PMC1239896](https://pmc.ncbi.nlm.nih.gov/articles/PMC1239896/)
  - JLinAlg methods: GSEA, gene-set enrichment.

<a id="wu-camera-2012"></a>
- **Di Wu and Gordon K. Smyth (2012).** Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Research 40:e133.* [DOI: 10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461) · [PMID: 22638577](https://pubmed.ncbi.nlm.nih.gov/22638577/) · [PMCID: PMC3458527](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/)
  - JLinAlg methods: CAMERA, correlation-aware enrichment.

<a id="deleeuw-magma-2015"></a>
- **Christiaan A. de Leeuw et al. (2015).** MAGMA: generalized gene-set analysis of GWAS data. *PLoS Computational Biology 11:e1004219.* [DOI: 10.1371/journal.pcbi.1004219](https://doi.org/10.1371/journal.pcbi.1004219) · [PMID: 25885710](https://pubmed.ncbi.nlm.nih.gov/25885710/) · [PMCID: PMC4401657](https://pmc.ncbi.nlm.nih.gov/articles/PMC4401657/)
  - JLinAlg methods: MAGMA, GWAS gene-set analysis.

## High-throughput inference

<a id="smyth-limma-2004"></a>
- **Gordon K. Smyth (2004).** Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology 3:Article 3.* [DOI: 10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027) · [PMID: 16646809](https://pubmed.ncbi.nlm.nih.gov/16646809/)
  - JLinAlg methods: moderated differential analysis, empirical-Bayes variance shrinkage.

<a id="law-voom-2014"></a>
- **Charity W. Law, Yunshun Chen, Wei Shi, and Gordon K. Smyth (2014).** voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology 15:R29.* [DOI: 10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29) · [PMID: 24485249](https://pubmed.ncbi.nlm.nih.gov/24485249/) · [PMCID: PMC4053721](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053721/)
  - JLinAlg methods: voom mean-variance precision weights, count differential analysis.

<a id="robinson-edger-2010"></a>
- **Mark D. Robinson, Davis J. McCarthy, and Gordon K. Smyth (2010).** edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics 26:139-140.* [DOI: 10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616) · [PMID: 19910308](https://pubmed.ncbi.nlm.nih.gov/19910308/) · [PMCID: PMC2796818](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/)
  - JLinAlg methods: negative-binomial differential analysis, dispersion moderation.

<a id="love-deseq2-2014"></a>
- **Michael I. Love, Wolfgang Huber, and Simon Anders (2014).** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology 15:550.* [DOI: 10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) · [PMID: 25516281](https://pubmed.ncbi.nlm.nih.gov/25516281/) · [PMCID: PMC4302049](https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/)
  - JLinAlg methods: median-ratio count normalization, dispersion shrinkage.

<a id="pedersen-combp-2012"></a>
- **Brent S. Pedersen, David A. Schwartz, Ivana V. Yang, and Katerina J. Kechris (2012).** Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values. *Bioinformatics 28:2986-2988.* [DOI: 10.1093/bioinformatics/bts545](https://doi.org/10.1093/bioinformatics/bts545) · [PMID: 22954632](https://pubmed.ncbi.nlm.nih.gov/22954632/) · [PMCID: PMC3496335](https://pmc.ncbi.nlm.nih.gov/articles/PMC3496335/)
  - JLinAlg methods: spatially correlated p-values, region-level genomic inference.

<a id="peters-dmrcate-2015"></a>
- **Timothy J. Peters et al. (2015).** De novo identification of differentially methylated regions in the human genome. *Epigenetics and Chromatin 8:6.* [DOI: 10.1186/1756-8935-8-6](https://doi.org/10.1186/1756-8935-8-6) · [PMID: 25972926](https://pubmed.ncbi.nlm.nih.gov/25972926/) · [PMCID: PMC4429355](https://pmc.ncbi.nlm.nih.gov/articles/PMC4429355/)
  - JLinAlg methods: differentially methylated regions, spatial methylation smoothing.

## Multiple testing

<a id="ignatiadis-ihw-2016"></a>
- **Nikolaos Ignatiadis, Bernd Klaus, Judith B. Zaugg, and Wolfgang Huber (2016).** Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nature Methods 13:577-580.* [DOI: 10.1038/nmeth.3885](https://doi.org/10.1038/nmeth.3885) · [PMID: 27240256](https://pubmed.ncbi.nlm.nih.gov/27240256/) · [PMCID: PMC4930141](https://pmc.ncbi.nlm.nih.gov/articles/PMC4930141/)
  - JLinAlg methods: independent hypothesis weighting, covariate-adaptive FDR.

<a id="meinshausen-hierarchy-2008"></a>
- **Nicolai Meinshausen (2008).** Hierarchical testing of variable importance. *Biometrika 95:265-278.* [DOI: 10.1093/biomet/asn007](https://doi.org/10.1093/biomet/asn007)
  - JLinAlg methods: hierarchical testing, familywise error control.

## Missing-data inference

<a id="rubin-mi-1987"></a>
- **Donald B. Rubin (1987).** Multiple Imputation for Nonresponse in Surveys. *Wiley.* [DOI: 10.1002/9780470316696](https://doi.org/10.1002/9780470316696)
  - JLinAlg methods: multiple-imputation pooling, missing-data uncertainty.

<a id="vanbuuren-mice-2011"></a>
- **Stef van Buuren and Karin Groothuis-Oudshoorn (2011).** mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software 45:1-67.* [DOI: 10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03)
  - JLinAlg methods: chained-equations imputation, mixed-type imputers.

<a id="barnard-rubin-1999"></a>
- **John Barnard and Donald B. Rubin (1999).** Small-sample degrees of freedom with multiple imputation. *Biometrika 86:948-955.* [DOI: 10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)
  - JLinAlg methods: Barnard-Rubin degrees of freedom, finite-sample multiple-imputation inference.

## Software interfaces

<a id="wilkinson-rogers-1973"></a>
- **G. N. Wilkinson and C. E. Rogers (1973).** Symbolic description of factorial models for analysis of variance. *Applied Statistics 22:392-399.* [DOI: 10.2307/2346786](https://doi.org/10.2307/2346786)
  - JLinAlg methods: statistical formulas, model matrices.

## Censored, ordinal, and sampling models

<a id="tobin-1958"></a>
- **James Tobin (1958).** Estimation of Relationships for Limited Dependent Variables. *Econometrica 26:24-36.* [DOI: 10.2307/1907382](https://doi.org/10.2307/1907382)
  - JLinAlg methods: Tobit, censored Gaussian regression.

<a id="kalbfleisch-prentice-2002"></a>
- **John D. Kalbfleisch and Ross L. Prentice (2002).** The Statistical Analysis of Failure Time Data, Second Edition. *Wiley Series in Probability and Statistics.* [DOI: 10.1002/9781118032985](https://doi.org/10.1002/9781118032985)
  - JLinAlg methods: parametric survival regression, accelerated failure time, censored likelihood.

<a id="mccullagh-1980"></a>
- **Peter McCullagh (1980).** Regression Models for Ordinal Data. *Journal of the Royal Statistical Society Series B 42:109-142.* [DOI: 10.1111/j.2517-6161.1980.tb01109.x](https://doi.org/10.1111/j.2517-6161.1980.tb01109.x)
  - JLinAlg methods: ordered logit, ordered probit, cumulative link models.

<a id="king-zeng-2001"></a>
- **Gary King and Langche Zeng (2001).** Logistic Regression in Rare Events Data. *Political Analysis 9:137-163.* [Primary source](https://gking.harvard.edu/files/0s.pdf)
  - JLinAlg methods: rare-events logistic regression, coefficient-bias correction.

<a id="king-zeng-isq-2001"></a>
- **Gary King and Langche Zeng (2001).** Explaining Rare Events in International Relations. *International Studies Quarterly 45:693-715.* [Primary source](https://gking.harvard.edu/files/gking/files/baby0s.pdf)
  - JLinAlg methods: case-control sampling, population-prevalence correction.

<a id="lumley-survey-2004"></a>
- **Thomas Lumley (2004).** Analysis of Complex Survey Samples. *Journal of Statistical Software 9(8):1-19.* [DOI: 10.18637/jss.v009.i08](https://doi.org/10.18637/jss.v009.i08)
  - JLinAlg methods: survey regression, sampling weights, stratified PSU Taylor covariance.

## Reproducible omics follow-up

<a id="sandve-reproducibility-2013"></a>
- **Geir Kjetil Sandve et al. (2013).** Ten Simple Rules for Reproducible Computational Research. *PLoS Computational Biology 9:e1003285.* [DOI: 10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285)
  - JLinAlg methods: Configuration provenance, Reproducible workflows.

<a id="wang-annovar-2010"></a>
- **Kai Wang, Mingyao Li, and Hakon Hakonarson (2010).** ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. *Nucleic Acids Research 38:e164.* [DOI: 10.1093/nar/gkq603](https://doi.org/10.1093/nar/gkq603)
  - JLinAlg methods: ANNOVAR consequence adapter.

<a id="mclaren-vep-2016"></a>
- **William McLaren et al. (2016).** The Ensembl Variant Effect Predictor. *Genome Biology 17:122.* [DOI: 10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)
  - JLinAlg methods: VEP consequence annotation, Variant consequence import.

<a id="langfelder-wgcna-2008"></a>
- **Peter Langfelder and Steve Horvath (2008).** WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics 9:559.* [DOI: 10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)
  - JLinAlg methods: Signed coexpression modules, Module eigengenes and hubs.

<a id="langfelder-preservation-2011"></a>
- **Peter Langfelder et al. (2011).** Is My Network Module Preserved and Reproducible?. *PLoS Computational Biology 7:e1001057.* [DOI: 10.1371/journal.pcbi.1001057](https://doi.org/10.1371/journal.pcbi.1001057)
  - JLinAlg methods: WGCNA module preservation.

<a id="meinshausen-networks-2006"></a>
- **Nicolai Meinshausen and Peter Buhlmann (2006).** High-dimensional graphs and variable selection with the Lasso. *Annals of Statistics 34:1436-1462.* [DOI: 10.1214/009053606000000281](https://doi.org/10.1214/009053606000000281)
  - JLinAlg methods: Gaussian neighborhood selection, Sparse conditional association.

<a id="huynh-thu-genie3-2010"></a>
- **Vân Anh Huynh-Thu et al. (2010).** Inferring Regulatory Networks from Expression Data Using Tree-Based Methods. *PLoS ONE 5:e12776.* [DOI: 10.1371/journal.pone.0012776](https://doi.org/10.1371/journal.pone.0012776)
  - JLinAlg methods: GENIE3 regulatory networks.

<a id="liu-carnival-2019"></a>
- **Anika Liu et al. (2019).** From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL. *npj Systems Biology and Applications 5:40.* [DOI: 10.1038/s41540-019-0118-z](https://doi.org/10.1038/s41540-019-0118-z)
  - JLinAlg methods: Inverse CARNIVAL, Contextual signaling networks.

<a id="fisher-correlation-1921"></a>
- **R. A. Fisher (1921).** On the probable error of a coefficient of correlation deduced from a small sample. *Metron 1:3-32.* [Primary source](https://digital.library.adelaide.edu.au/dspace/handle/2440/15169)
  - JLinAlg methods: Fisher correlation transformation, Independent-group differential correlation.

<a id="brin-pagerank-1998"></a>
- **Sergey Brin and Lawrence Page (1998).** The anatomy of a large-scale hypertextual Web search engine. *Computer Networks and ISDN Systems 30:107-117.* [Primary source](https://snap.stanford.edu/class/cs224w-readings/Brin98Anatomy.pdf)
  - JLinAlg methods: Random walk with restart, Network propagation.
