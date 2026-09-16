# Scientific citations

This bibliography links JLinAlg features to the primary publications that introduced or established their statistical methods. Cite both the relevant method paper and JLinAlg when reporting an analysis. Inclusion here documents methodological provenance; it is not a claim that every implementation detail is identical to the cited software.

The website provides a [searchable citation index](https://robbyjo.github.io/JLinAlg/citations.html). This file and the per-vignette citation panels are generated from [citations.json](citations.json).

## Foundational inference

<a id="student-1908"></a>
### Student [W. S. Gosset] (1908)

**The probable error of a mean.** Biometrika 6:1-25. [Primary source](https://doi.org/10.2307/2331554)

JLinAlg features: Student t test, small-sample mean inference.

<a id="welch-1947"></a>
### B. L. Welch (1947)

**The generalization of Student's problem when several different population variances are involved.** Biometrika 34:28-35. [Primary source](https://doi.org/10.1093/biomet/34.1-2.28)

JLinAlg features: Welch t test, unequal-variance inference.

<a id="pearson-1895"></a>
### Karl Pearson (1895)

**Notes on regression and inheritance in the case of two parents.** Proceedings of the Royal Society of London 58:240-242. [Primary source](https://doi.org/10.1098/rspl.1895.0041)

JLinAlg features: Pearson correlation, linear association.

<a id="spearman-1904"></a>
### Charles Spearman (1904)

**The proof and measurement of association between two things.** American Journal of Psychology 15:72-101. [Primary source](https://doi.org/10.2307/1412159)

JLinAlg features: Spearman correlation, rank association.

<a id="wilcoxon-1945"></a>
### Frank Wilcoxon (1945)

**Individual comparisons by ranking methods.** Biometrics Bulletin 1:80-83. [Primary source](https://doi.org/10.2307/3001968)

JLinAlg features: Wilcoxon signed-rank test, rank inference.

<a id="mann-whitney-1947"></a>
### H. B. Mann and D. R. Whitney (1947)

**On a test of whether one of two random variables is stochastically larger than the other.** Annals of Mathematical Statistics 18:50-60. [Primary source](https://doi.org/10.1214/aoms/1177730491)

JLinAlg features: Mann-Whitney U test, two-sample rank inference.

<a id="fisher-1922"></a>
### R. A. Fisher (1922)

**On the interpretation of chi-square from contingency tables, and the calculation of P.** Journal of the Royal Statistical Society 85:87-94. [Primary source](https://doi.org/10.2307/2340521)

JLinAlg features: Fisher exact test, contingency tables.

<a id="mcnemar-1947"></a>
### Quinn McNemar (1947)

**Note on the sampling error of the difference between correlated proportions or percentages.** Psychometrika 12:153-157. [Primary source](https://doi.org/10.1007/BF02295996)

JLinAlg features: McNemar test, paired binary data.

<a id="friedman-1937"></a>
### Milton Friedman (1937)

**The use of ranks to avoid the assumption of normality implicit in the analysis of variance.** Journal of the American Statistical Association 32:675-701. [Primary source](https://doi.org/10.1080/01621459.1937.10503522)

JLinAlg features: Friedman test, blocked rank inference.

<a id="benjamini-hochberg-1995"></a>
### Yoav Benjamini and Yosef Hochberg (1995)

**Controlling the false discovery rate: a practical and powerful approach to multiple testing.** Journal of the Royal Statistical Society B 57:289-300. [Primary source](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

JLinAlg features: Benjamini-Hochberg FDR, multiple testing.

## Regression

<a id="nelder-wedderburn-1972"></a>
### J. A. Nelder and R. W. M. Wedderburn (1972)

**Generalized linear models.** Journal of the Royal Statistical Society A 135:370-384. [Primary source](https://doi.org/10.2307/2344614)

JLinAlg features: generalized linear models, IRLS, logistic regression, Poisson regression.

<a id="bliss-1934"></a>
### C. I. Bliss (1934)

**The method of probits.** Science 79:38-39. [Primary source](https://doi.org/10.1126/science.79.2037.38)

JLinAlg features: probit regression, binary response.

<a id="hoerl-kennard-1970"></a>
### Arthur E. Hoerl and Robert W. Kennard (1970)

**Ridge regression: biased estimation for nonorthogonal problems.** Technometrics 12:55-67. [Primary source](https://doi.org/10.1080/00401706.1970.10488634)

JLinAlg features: ridge regression, prediction scores.

<a id="tibshirani-1996"></a>
### Robert Tibshirani (1996)

**Regression shrinkage and selection via the lasso.** Journal of the Royal Statistical Society B 58:267-288. [Primary source](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x)

JLinAlg features: lasso, sparse regression.

<a id="zou-hastie-2005"></a>
### Hui Zou and Trevor Hastie (2005)

**Regularization and variable selection via the elastic net.** Journal of the Royal Statistical Society B 67:301-320. [Primary source](https://doi.org/10.1111/j.1467-9868.2005.00503.x)

JLinAlg features: elastic net, grouped shrinkage.

<a id="lee-2016"></a>
### Jason D. Lee, Dennis L. Sun, Yuekai Sun, and Jonathan E. Taylor (2016)

**Exact post-selection inference, with application to the lasso.** Annals of Statistics 44:907-927. [Primary source](https://doi.org/10.1214/15-AOS1371)

JLinAlg features: polyhedral selective inference, lasso inference.

<a id="koenker-bassett-1978"></a>
### Roger Koenker and Gilbert Bassett Jr. (1978)

**Regression quantiles.** Econometrica 46:33-50. [Primary source](https://doi.org/10.2307/1913643)

JLinAlg features: quantile regression.

<a id="friedman-1984"></a>
### Jerome H. Friedman (1984)

**A variable span smoother.** Laboratory for Computational Statistics Technical Report 5. [Primary source](https://doi.org/10.2172/1447470)

JLinAlg features: supersmoother, adaptive span smoothing.

<a id="nadaraya-1964"></a>
### E. A. Nadaraya (1964)

**On estimating regression.** Theory of Probability and Its Applications 9:141-142. [Primary source](https://doi.org/10.1137/1109020)

JLinAlg features: kernel regression, Nadaraya-Watson estimator.

<a id="robinson-1988"></a>
### Peter M. Robinson (1988)

**Root-N-consistent semiparametric regression.** Econometrica 56:931-954. [Primary source](https://doi.org/10.2307/1912705)

JLinAlg features: partially linear regression, semiparametric regression.

<a id="cleveland-1979"></a>
### William S. Cleveland (1979)

**Robust locally weighted regression and smoothing scatterplots.** Journal of the American Statistical Association 74:829-836. [Primary source](https://doi.org/10.1080/01621459.1979.10481038)

JLinAlg features: LOESS, robust local regression.

<a id="hastie-tibshirani-1986"></a>
### Trevor Hastie and Robert Tibshirani (1986)

**Generalized additive models.** Statistical Science 1:297-310. [Primary source](https://doi.org/10.1214/ss/1177013604)

JLinAlg features: generalized additive models, additive smooths.

<a id="eilers-marx-1996"></a>
### Paul H. C. Eilers and Brian D. Marx (1996)

**Flexible smoothing with B-splines and penalties.** Statistical Science 11:89-121. [Primary source](https://doi.org/10.1214/ss/1038425655)

JLinAlg features: P-splines, penalized B-splines.

<a id="wood-2011"></a>
### Simon N. Wood (2011)

**Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models.** Journal of the Royal Statistical Society B 73:3-36. [Primary source](https://doi.org/10.1111/j.1467-9868.2010.00749.x)

JLinAlg features: GAM REML, smoothness selection, penalized GLM.

<a id="rigby-stasinopoulos-2005"></a>
### Robert A. Rigby and D. Mikis Stasinopoulos (2005)

**Generalized additive models for location, scale and shape.** Applied Statistics 54:507-554. [Primary source](https://doi.org/10.1111/j.1467-9876.2005.00510.x)

JLinAlg features: GAMLSS, distributional regression.

<a id="yee-2010"></a>
### Thomas W. Yee (2010)

**The VGAM package for categorical data analysis.** Journal of Statistical Software 32:1-34. [Primary source](https://doi.org/10.18637/jss.v032.i10)

JLinAlg features: VGAM, ordinal and categorical regression.

<a id="ferrari-cribari-neto-2004"></a>
### Silvia Ferrari and Francisco Cribari-Neto (2004)

**Beta regression for modelling rates and proportions.** Journal of Applied Statistics 31:799-815. [Primary source](https://doi.org/10.1080/0266476042000214501)

JLinAlg features: beta regression, mean-precision models.

<a id="beasley-2009"></a>
### T. Mark Beasley, Stephen Erickson, and David B. Allison (2009)

**Rank-based inverse normal transformations are increasingly used, but are they merited?.** Behavior Genetics 39:580-595. [Primary source](https://doi.org/10.1007/s10519-009-9281-0)

JLinAlg features: rank inverse-normal transform, omics transforms.

## Mixed, longitudinal, and survival models

<a id="patterson-thompson-1971"></a>
### H. D. Patterson and Robin Thompson (1971)

**Recovery of inter-block information when block sizes are unequal.** Biometrika 58:545-554. [Primary source](https://doi.org/10.1093/biomet/58.3.545)

JLinAlg features: REML, variance components.

<a id="henderson-1975"></a>
### C. R. Henderson (1975)

**Best linear unbiased estimation and prediction under a selection model.** Biometrics 31:423-447. [Primary source](https://doi.org/10.2307/2529430)

JLinAlg features: BLUP, mixed-model equations.

<a id="kenward-roger-1997"></a>
### Michael G. Kenward and James H. Roger (1997)

**Small sample inference for fixed effects from restricted maximum likelihood.** Biometrics 53:983-997. [Primary source](https://doi.org/10.2307/2533558)

JLinAlg features: Kenward-Roger inference, small-sample mixed models.

<a id="breslow-clayton-1993"></a>
### Norman E. Breslow and David G. Clayton (1993)

**Approximate inference in generalized linear mixed models.** Journal of the American Statistical Association 88:9-25. [Primary source](https://doi.org/10.1080/01621459.1993.10594284)

JLinAlg features: PQL, generalized linear mixed models.

<a id="lambert-1992"></a>
### Diane Lambert (1992)

**Zero-inflated Poisson regression, with an application to defects in manufacturing.** Technometrics 34:1-14. [Primary source](https://doi.org/10.1080/00401706.1992.10485228)

JLinAlg features: zero-inflated Poisson, mixture count models.

<a id="lindstrom-bates-1990"></a>
### Mary J. Lindstrom and Douglas M. Bates (1990)

**Nonlinear mixed effects models for repeated measures data.** Biometrics 46:673-687. [Primary source](https://doi.org/10.2307/2532087)

JLinAlg features: nonlinear mixed-effects models.

<a id="liang-zeger-1986"></a>
### Kung-Yee Liang and Scott L. Zeger (1986)

**Longitudinal data analysis using generalized linear models.** Biometrika 73:13-22. [Primary source](https://doi.org/10.1093/biomet/73.1.13)

JLinAlg features: GEE, working correlation, sandwich inference.

<a id="cox-1972"></a>
### D. R. Cox (1972)

**Regression models and life-tables.** Journal of the Royal Statistical Society B 34:187-220. [Primary source](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)

JLinAlg features: Cox proportional hazards, partial likelihood.

<a id="efron-1977"></a>
### Bradley Efron (1977)

**The efficiency of Cox's likelihood function for censored data.** Journal of the American Statistical Association 72:557-565. [Primary source](https://doi.org/10.1080/01621459.1977.10480613)

JLinAlg features: Efron ties, Cox regression.

## Quantitative and statistical genetics

<a id="henderson-1976"></a>
### C. R. Henderson (1976)

**A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values.** Biometrics 32:69-83. [Primary source](https://doi.org/10.2307/2529339)

JLinAlg features: pedigree A inverse, animal models.

<a id="vanraden-2008"></a>
### Paul M. VanRaden (2008)

**Efficient methods to compute genomic predictions.** Journal of Dairy Science 91:4414-4423. [Primary source](https://doi.org/10.3168/jds.2007-0980)

JLinAlg features: genomic relationship matrix, genomic prediction.

<a id="yang-gcta-2011"></a>
### Jian Yang et al. (2011)

**GCTA: a tool for genome-wide complex trait analysis.** American Journal of Human Genetics 88:76-82. [Primary source](https://doi.org/10.1016/j.ajhg.2010.11.011)

JLinAlg features: GRM, SNP heritability, GCTA.

<a id="kang-emmax-2010"></a>
### Hyun Min Kang et al. (2010)

**Variance component model to account for sample structure in genome-wide association studies.** Nature Genetics 42:348-354. [Primary source](https://doi.org/10.1038/ng.548)

JLinAlg features: EMMAX, mixed-model GWAS.

<a id="zhang-p3d-2010"></a>
### Zhiwu Zhang et al. (2010)

**Mixed linear model approach adapted for genome-wide association studies.** Nature Genetics 42:355-360. [Primary source](https://doi.org/10.1038/ng.546)

JLinAlg features: P3D, compressed mixed-model GWAS.

<a id="wu-skat-2011"></a>
### Michael C. Wu et al. (2011)

**Rare-variant association testing for sequencing data with the sequence kernel association test.** American Journal of Human Genetics 89:82-93. [Primary source](https://doi.org/10.1016/j.ajhg.2011.05.029)

JLinAlg features: SKAT, rare-variant set tests.

<a id="lee-skato-2012"></a>
### Seunggeun Lee et al. (2012)

**Optimal unified approach for rare-variant association testing.** American Journal of Human Genetics 91:224-237. [Primary source](https://doi.org/10.1016/j.ajhg.2012.06.007)

JLinAlg features: SKAT-O, burden-kernel omnibus tests.

<a id="madsen-browning-2009"></a>
### B. E. Madsen and S. R. Browning (2009)

**A groupwise association test for rare mutations using a weighted sum statistic.** PLoS Genetics 5:e1000384. [Primary source](https://doi.org/10.1371/journal.pgen.1000384)

JLinAlg features: weighted burden test, rare variants.

<a id="yang-cojo-2012"></a>
### Jian Yang et al. (2012)

**Conditional and joint multiple-SNP analysis of GWAS summary statistics.** Nature Genetics 44:369-375. [Primary source](https://doi.org/10.1038/ng.2213)

JLinAlg features: COJO, conditional summary statistics.

<a id="bulik-sullivan-ldsc-2015"></a>
### Brendan K. Bulik-Sullivan et al. (2015)

**LD Score regression distinguishes confounding from polygenicity in genome-wide association studies.** Nature Genetics 47:291-295. [Primary source](https://doi.org/10.1038/ng.3211)

JLinAlg features: LD Score regression, SNP heritability, confounding intercept.

<a id="bulik-sullivan-rg-2015"></a>
### Brendan K. Bulik-Sullivan et al. (2015)

**An atlas of genetic correlations across human diseases and traits.** Nature Genetics 47:1236-1241. [Primary source](https://doi.org/10.1038/ng.3406)

JLinAlg features: genetic correlation, bivariate LDSC.

<a id="gamazon-predixcan-2015"></a>
### Eric R. Gamazon et al. (2015)

**A gene-based association method for mapping traits using reference transcriptome data.** Nature Genetics 47:1091-1098. [Primary source](https://doi.org/10.1038/ng.3367)

JLinAlg features: PrediXcan, genetically predicted expression, TWAS.

<a id="gusev-fusion-2016"></a>
### Alexander Gusev et al. (2016)

**Integrative approaches for large-scale transcriptome-wide association studies.** Nature Genetics 48:245-252. [Primary source](https://doi.org/10.1038/ng.3506)

JLinAlg features: FUSION TWAS, summary TWAS, PWAS.

<a id="grotzinger-genomic-sem-2019"></a>
### Andrew D. Grotzinger et al. (2019)

**Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits.** Nature Human Behaviour 3:513-525. [Primary source](https://doi.org/10.1038/s41562-019-0566-x)

JLinAlg features: Genomic SEM, shared genetic factors.

<a id="euesden-prsice-2015"></a>
### Jack Euesden, Cathryn M. Lewis, and Paul F. O'Reilly (2015)

**PRSice: polygenic risk score software.** Bioinformatics 31:1466-1468. [Primary source](https://doi.org/10.1093/bioinformatics/btu848)

JLinAlg features: polygenic scores, score evaluation.

## Latent variation and batch correction

<a id="leek-storey-2007"></a>
### Jeffrey T. Leek and John D. Storey (2007)

**Capturing heterogeneity in gene expression studies by surrogate variable analysis.** PLoS Genetics 3:e161. [Primary source](https://doi.org/10.1371/journal.pgen.0030161)

JLinAlg features: SVA, latent confounders.

<a id="stegle-peer-2012"></a>
### Oliver Stegle et al. (2012)

**Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses.** Nature Protocols 7:500-507. [Primary source](https://doi.org/10.1038/nprot.2011.457)

JLinAlg features: PEER, Bayesian factor analysis.

<a id="johnson-combat-2007"></a>
### W. Evan Johnson, Cheng Li, and Ariel Rabinovic (2007)

**Adjusting batch effects in microarray expression data using empirical Bayes methods.** Biostatistics 8:118-127. [Primary source](https://doi.org/10.1093/biostatistics/kxj037)

JLinAlg features: ComBat, batch correction, empirical Bayes.

## Causal and mediation analysis

<a id="baron-kenny-1986"></a>
### Reuben M. Baron and David A. Kenny (1986)

**The moderator-mediator variable distinction in social psychological research.** Journal of Personality and Social Psychology 51:1173-1182. [Primary source](https://doi.org/10.1037/0022-3514.51.6.1173)

JLinAlg features: mediation, path decomposition.

<a id="sobel-1982"></a>
### Michael E. Sobel (1982)

**Asymptotic confidence intervals for indirect effects in structural equation models.** Sociological Methodology 13:290-312. [Primary source](https://doi.org/10.2307/270723)

JLinAlg features: Sobel test, indirect-effect inference.

<a id="staiger-stock-1997"></a>
### Douglas Staiger and James H. Stock (1997)

**Instrumental variables regression with weak instruments.** Econometrica 65:557-586. [Primary source](https://doi.org/10.2307/2171753)

JLinAlg features: 2SLS, weak instruments, first-stage diagnostics.

<a id="white-1980"></a>
### Halbert White (1980)

**A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity.** Econometrica 48:817-838. [Primary source](https://doi.org/10.2307/1912934)

JLinAlg features: heteroskedasticity-robust covariance, robust 2SLS inference.

<a id="bowden-mr-egger-2015"></a>
### Jack Bowden, George Davey Smith, and Stephen Burgess (2015)

**Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression.** International Journal of Epidemiology 44:512-525. [Primary source](https://doi.org/10.1093/ije/dyv080)

JLinAlg features: MR-Egger, directional pleiotropy, Egger intercept.

<a id="bowden-weighted-median-2016"></a>
### Jack Bowden et al. (2016)

**Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator.** Genetic Epidemiology 40:304-314. [Primary source](https://doi.org/10.1002/gepi.21965)

JLinAlg features: weighted-median MR, invalid instruments.

<a id="hemani-steiger-2017"></a>
### Gibran Hemani, Kate Tilling, and George Davey Smith (2017)

**Orienting the causal relationship between imprecisely measured traits using GWAS summary data.** PLoS Genetics 13:e1007081. [Primary source](https://doi.org/10.1371/journal.pgen.1007081)

JLinAlg features: Steiger directionality, MR orientation.

<a id="zhao-mr-raps-2020"></a>
### Qingyuan Zhao et al. (2020)

**Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score.** Annals of Statistics 48:1742-1769. [Primary source](https://doi.org/10.1214/19-AOS1866)

JLinAlg features: MR-RAPS, weak instruments, idiosyncratic pleiotropy.

<a id="verbanck-mr-presso-2018"></a>
### Marie Verbanck et al. (2018)

**Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization.** Nature Genetics 50:693-698. [Primary source](https://doi.org/10.1038/s41588-018-0099-7)

JLinAlg features: MR-PRESSO, pleiotropic outliers.

<a id="burgess-mvmr-2015"></a>
### Stephen Burgess and Simon G. Thompson (2015)

**Multivariable Mendelian randomization: the use of pleiotropic genetic variants to estimate causal effects.** American Journal of Epidemiology 181:251-260. [Primary source](https://doi.org/10.1093/aje/kwu283)

JLinAlg features: multivariable MR, direct causal effects.

<a id="burgess-overlap-2016"></a>
### Stephen Burgess, Neil M. Davies, and Simon G. Thompson (2016)

**Bias due to participant overlap in two-sample Mendelian randomization.** Genetic Epidemiology 40:597-608. [Primary source](https://doi.org/10.1002/gepi.21998)

JLinAlg features: sample overlap, two-sample MR.

<a id="hemani-mrbase-2018"></a>
### Gibran Hemani et al. (2018)

**The MR-Base platform supports systematic causal inference across the human phenome.** eLife 7:e34408. [Primary source](https://doi.org/10.7554/eLife.34408)

JLinAlg features: TwoSampleMR, MR-Base, phenome-wide MR.

## Evidence synthesis

<a id="dersimonian-laird-1986"></a>
### Rebecca DerSimonian and Nan Laird (1986)

**Meta-analysis in clinical trials.** Controlled Clinical Trials 7:177-188. [Primary source](https://doi.org/10.1016/0197-2456(86)90046-2)

JLinAlg features: random-effects meta-analysis, DerSimonian-Laird.

<a id="paule-mandel-1982"></a>
### Robert C. Paule and John Mandel (1982)

**Consensus values and weighting factors.** Journal of Research of the National Bureau of Standards 87:377-385. [Primary source](https://doi.org/10.6028/jres.087.022)

JLinAlg features: Paule-Mandel tau squared, random-effects meta-analysis.

<a id="knapp-hartung-2003"></a>
### Guido Knapp and Joachim Hartung (2003)

**Improved tests for a random effects meta-regression with a single covariate.** Statistics in Medicine 22:2693-2710. [Primary source](https://doi.org/10.1002/sim.1482)

JLinAlg features: Knapp-Hartung inference, meta-regression.

<a id="higgins-thompson-2002"></a>
### Julian P. T. Higgins and Simon G. Thompson (2002)

**Quantifying heterogeneity in a meta-analysis.** Statistics in Medicine 21:1539-1558. [Primary source](https://doi.org/10.1002/sim.1186)

JLinAlg features: I-squared, heterogeneity.

## Time series

<a id="kalman-1960"></a>
### Rudolf E. Kalman (1960)

**A new approach to linear filtering and prediction problems.** Journal of Basic Engineering 82:35-45. [Primary source](https://doi.org/10.1115/1.3662552)

JLinAlg features: Kalman filter, state-space models, diffuse smoothing.

<a id="ljung-box-1978"></a>
### Greta M. Ljung and George E. P. Box (1978)

**On a measure of lack of fit in time series models.** Biometrika 65:297-303. [Primary source](https://doi.org/10.1093/biomet/65.2.297)

JLinAlg features: Ljung-Box test, residual diagnostics.

<a id="hyndman-khandakar-2008"></a>
### Rob J. Hyndman and Yeasmin Khandakar (2008)

**Automatic time series forecasting: the forecast package for R.** Journal of Statistical Software 27:1-22. [Primary source](https://doi.org/10.18637/jss.v027.i03)

JLinAlg features: automatic ARIMA, stepwise order selection.

## Fine mapping, colocalization, and SEM

<a id="wang-susie-2020"></a>
### Gao Wang, Abhishek Sarkar, Peter Carbonetto, and Matthew Stephens (2020)

**A simple new approach to variable selection in regression, with application to genetic fine mapping.** Journal of the Royal Statistical Society B 82:1273-1300. [Primary source](https://doi.org/10.1111/rssb.12388)

JLinAlg features: SuSiE, credible sets, fine mapping.

<a id="giambartolomei-coloc-2014"></a>
### Claudia Giambartolomei et al. (2014)

**Bayesian test for colocalisation between pairs of genetic association studies using summary statistics.** PLoS Genetics 10:e1004383. [Primary source](https://doi.org/10.1371/journal.pgen.1004383)

JLinAlg features: coloc, H0-H4 hypotheses.

<a id="wallace-coloc-2021"></a>
### Chris Wallace (2021)

**A more accurate method for colocalisation analysis allowing for multiple causal variants.** PLoS Genetics 17:e1009440. [Primary source](https://doi.org/10.1371/journal.pgen.1009440)

JLinAlg features: SuSiE colocalization, multiple causal variants.

<a id="wright-1921"></a>
### Sewall Wright (1921)

**Correlation and causation.** Journal of Agricultural Research 20:557-585. [Primary source](https://naldc.nal.usda.gov/catalog/IND43966364)

JLinAlg features: path analysis, structural equation models.

<a id="mcardle-mcdonald-1984"></a>
### John J. McArdle and Roderick P. McDonald (1984)

**Some algebraic properties of the Reticular Action Model for moment structures.** British Journal of Mathematical and Statistical Psychology 37:234-251. [Primary source](https://doi.org/10.1111/j.2044-8317.1984.tb00802.x)

JLinAlg features: RAM formulation, SEM covariance structure.

<a id="enders-bandalos-2001"></a>
### Craig K. Enders and Deborah L. Bandalos (2001)

**The relative performance of full information maximum likelihood estimation for missing data in structural equation models.** Structural Equation Modeling 8:430-457. [Primary source](https://doi.org/10.1207/S15328007SEM0803_5)

JLinAlg features: FIML, SEM missing data.

<a id="katsikatsou-2012"></a>
### Myrsini Katsikatsou et al. (2012)

**Pairwise likelihood estimation for factor analysis models with ordinal data.** Computational Statistics and Data Analysis 56:4243-4258. [Primary source](https://doi.org/10.1016/j.csda.2012.04.010)

JLinAlg features: ordinal pairwise likelihood, ordinal SEM.

## Enrichment and annotation

<a id="subramanian-gsea-2005"></a>
### Aravind Subramanian et al. (2005)

**Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.** Proceedings of the National Academy of Sciences 102:15545-15550. [Primary source](https://doi.org/10.1073/pnas.0506580102)

JLinAlg features: GSEA, gene-set enrichment.

<a id="wu-camera-2012"></a>
### Di Wu and Gordon K. Smyth (2012)

**Camera: a competitive gene set test accounting for inter-gene correlation.** Nucleic Acids Research 40:e133. [Primary source](https://doi.org/10.1093/nar/gks461)

JLinAlg features: CAMERA, correlation-aware enrichment.

<a id="deleeuw-magma-2015"></a>
### Christiaan A. de Leeuw et al. (2015)

**MAGMA: generalized gene-set analysis of GWAS data.** PLoS Computational Biology 11:e1004219. [Primary source](https://doi.org/10.1371/journal.pcbi.1004219)

JLinAlg features: MAGMA, GWAS gene-set analysis.

## Software interfaces

<a id="wilkinson-rogers-1973"></a>
### G. N. Wilkinson and C. E. Rogers (1973)

**Symbolic description of factorial models for analysis of variance.** Applied Statistics 22:392-399. [Primary source](https://doi.org/10.2307/2346786)

JLinAlg features: statistical formulas, model matrices.
