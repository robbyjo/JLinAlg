# JLinAlg cross-language numerical reference generator.
#
# The committed Java tests contain the printed values from this script. R is
# not required to run the Java test suite. To regenerate, install R plus:
#   nlme 3.1-169, lme4, pbkrtest, rrBLUP 4.6.3,
#   MendelianRandomization 0.10.0

options(digits = 17, contrasts = c("contr.treatment", "contr.poly"))
arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
script_directory <- dirname(normalizePath(script_argument))
project_root <- normalizePath(file.path(script_directory, "../../../.."))
reference_library <- Sys.getenv(
  "JLINALG_R_REFERENCE_LIBRARY",
  file.path(project_root, ".r-reference-lib")
)
.libPaths(c(reference_library, .libPaths()))

cat("R.version=", R.version.string, "\n", sep = "")
cat("stats.version=", as.character(packageVersion("stats")), "\n", sep = "")
cat("nlme.version=", as.character(packageVersion("nlme")), "\n", sep = "")
cat("lme4.version=", as.character(packageVersion("lme4")), "\n", sep = "")
cat("pbkrtest.version=", as.character(packageVersion("pbkrtest")), "\n", sep = "")
cat("rrBLUP.version=", as.character(packageVersion("rrBLUP")), "\n", sep = "")
cat("MendelianRandomization.version=",
    as.character(packageVersion("MendelianRandomization")), "\n", sep = "")

emit <- function(name, value) {
  cat(name, "=")
  dput(unname(value))
}

# Base R documented data: ordinary Gaussian regression.
lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
lm_summary <- summary(lm_fit)
emit("lm.mtcars.coefficients", coef(lm_fit))
emit("lm.mtcars.standard_errors", lm_summary$coefficients[, "Std. Error"])
emit("lm.mtcars.t_statistics", lm_summary$coefficients[, "t value"])
emit("lm.mtcars.p_values", lm_summary$coefficients[, "Pr(>|t|)"])
emit("lm.mtcars.rss", deviance(lm_fit))
emit("lm.mtcars.residual_variance", lm_summary$sigma^2)
emit("lm.mtcars.log_likelihood", as.numeric(logLik(lm_fit)))
emit("mtcars.mpg", mtcars$mpg)
emit("mtcars.wt", mtcars$wt)
emit("mtcars.hp", mtcars$hp)
emit("mtcars.am", mtcars$am)

# Base R binomial GLM on the same documented data.
binomial_fit <- glm(am ~ wt + hp, data = mtcars, family = binomial())
binomial_summary <- summary(binomial_fit)
emit("glm.mtcars.binomial.coefficients", coef(binomial_fit))
emit("glm.mtcars.binomial.standard_errors",
     binomial_summary$coefficients[, "Std. Error"])
emit("glm.mtcars.binomial.deviance", deviance(binomial_fit))
emit("glm.mtcars.binomial.log_likelihood", as.numeric(logLik(binomial_fit)))
emit("glm.mtcars.binomial.aic", AIC(binomial_fit))

# Base R's documented warpbreaks data and a canonical Poisson GLM.
poisson_fit <- glm(breaks ~ wool + tension,
                   data = warpbreaks, family = poisson())
poisson_summary <- summary(poisson_fit)
emit("glm.warpbreaks.poisson.coefficients", coef(poisson_fit))
emit("glm.warpbreaks.poisson.standard_errors",
     poisson_summary$coefficients[, "Std. Error"])
emit("glm.warpbreaks.poisson.deviance", deviance(poisson_fit))
emit("glm.warpbreaks.poisson.log_likelihood", as.numeric(logLik(poisson_fit)))
emit("glm.warpbreaks.poisson.aic", AIC(poisson_fit))
emit("warpbreaks.breaks", warpbreaks$breaks)

# Base R conditional AR, ARMA, and seasonal ARIMA reference fits.
ar_lh <- arima(lh, order = c(2, 0, 0), method = "CSS")
arma_lh <- arima(lh, order = c(1, 0, 1), method = "CSS")
seasonal_deaths <- arima(USAccDeaths, order = c(0, 1, 1),
                         seasonal = list(order = c(0, 1, 1), period = 12),
                         method = "CSS")
emit("stats.lh.ar2.coefficients", coef(ar_lh))
emit("stats.lh.ar2.innovation_variance", ar_lh$sigma2)
emit("stats.lh.arma11.coefficients", coef(arma_lh))
emit("stats.lh.arma11.innovation_variance", arma_lh$sigma2)
emit("stats.USAccDeaths.sarima.coefficients", coef(seasonal_deaths))
emit("stats.USAccDeaths.sarima.innovation_variance", seasonal_deaths$sigma2)

# nlme manual example: random-intercept REML on ergoStool.
suppressPackageStartupMessages(library(nlme))
ergo <- as.data.frame(ergoStool)
ergo$Type <- factor(as.character(ergo$Type), levels = c("T1", "T2", "T3", "T4"))
nlme_fit <- lme(effort ~ Type, data = ergo,
                random = ~ 1 | Subject, method = "REML")
emit("nlme.ergoStool.fixed_effects", fixef(nlme_fit))
emit("nlme.ergoStool.fixed_covariance", as.vector(t(vcov(nlme_fit))))
emit("nlme.ergoStool.random_variance",
     as.numeric(VarCorr(nlme_fit)[1, "Variance"]))
emit("nlme.ergoStool.residual_variance", sigma(nlme_fit)^2)
emit("nlme.ergoStool.log_likelihood", as.numeric(logLik(nlme_fit)))
emit("nlme.ergoStool.random_effects", ranef(nlme_fit)[, 1])
emit("ergoStool.effort", ergo$effort)
emit("ergoStool.subject", as.integer(as.character(ergo$Subject)))
emit("ergoStool.type", as.integer(ergo$Type))

# lme4/pbkrtest Kenward-Roger reference on the same random-intercept model.
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(pbkrtest))
lme4_fit <- lmer(effort ~ Type + (1 | Subject), data = ergo, REML = TRUE)
kr_unadjusted <- vcov(lme4_fit)
kr_adjusted <- vcovAdj(lme4_fit)
kr_identity <- diag(ncol(model.matrix(lme4_fit)))
kr_df <- vapply(seq_len(nrow(kr_identity)), function(index) {
  Lb_ddf(kr_identity[index, ], kr_unadjusted, kr_adjusted)
}, numeric(1))
emit("pbkrtest.ergoStool.fixed_covariance_adjusted",
     as.vector(t(kr_adjusted)))
emit("pbkrtest.ergoStool.coefficient_df", kr_df)

# rrBLUP reference for a pedigree-like additive relationship matrix with
# unphenotyped founders and repeated records.
suppressPackageStartupMessages(library(rrBLUP))
parents <- matrix(c(
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  1, 2,
  1, 2,
  3, 4,
  3, 4,
  5, 7,
  6, 8
), ncol = 2, byrow = TRUE)
A <- matrix(0, 10, 10)
for (i in seq_len(10)) {
  sire <- parents[i, 1]
  dam <- parents[i, 2]
  if (i > 1) {
    for (j in seq_len(i - 1)) {
      A[i, j] <- A[j, i] <- 0.5 *
        ((if (sire == 0) 0 else A[sire, j]) +
         (if (dam == 0) 0 else A[dam, j]))
    }
  }
  A[i, i] <- 1 + if (sire == 0 || dam == 0) 0 else 0.5 * A[sire, dam]
}
observed_animals <- rep(5:10, each = 2)
pedigree_y <- c(8.5, 9.1, 8.9, 9.5, 10.8, 11.4,
                10.4, 11.0, 9.3, 9.9, 10.0, 10.6)
Z <- matrix(0, length(pedigree_y), 10)
Z[cbind(seq_along(pedigree_y), observed_animals)] <- 1
rrblup_fit <- mixed.solve(y = pedigree_y, Z = Z, K = A,
                          X = matrix(1, length(pedigree_y), 1),
                          method = "REML", SE = TRUE)
emit("rrBLUP.pedigree.genetic_variance", rrblup_fit$Vu)
emit("rrBLUP.pedigree.residual_variance", rrblup_fit$Ve)
emit("rrBLUP.pedigree.fixed_effect", rrblup_fit$beta)
emit("rrBLUP.pedigree.breeding_values", rrblup_fit$u)
emit("rrBLUP.pedigree.prediction_error_variances", rrblup_fit$u.SE^2)
emit("rrBLUP.pedigree.log_likelihood", rrblup_fit$LL)
emit("rrBLUP.pedigree.relationship", as.vector(t(A)))

# CRAN MendelianRandomization manual example using its shipped LDL-C data.
suppressPackageStartupMessages(library(MendelianRandomization))
mr_data <- mr_input(bx = ldlc, bxse = ldlcse,
                    by = chdlodds, byse = chdloddsse)
mr_ivw_fixed <- mr_ivw(mr_data, model = "fixed",
                       weights = "simple", distribution = "normal")
mr_ivw_random <- mr_ivw(mr_data, model = "random",
                        weights = "simple", distribution = "normal")
mr_egger_random <- mr_egger(mr_data, distribution = "normal")
set.seed(20260831)
mr_weighted_median <- mr_median(mr_data, weighting = "weighted",
                                iterations = 10000)
emit("mr.ldlc.bx", ldlc)
emit("mr.ldlc.bxse", ldlcse)
emit("mr.ldlc.by", chdlodds)
emit("mr.ldlc.byse", chdloddsse)
emit("mr.ldlc.ivw_fixed.estimate", mr_ivw_fixed@Estimate)
emit("mr.ldlc.ivw_fixed.standard_error", mr_ivw_fixed@StdError)
emit("mr.ldlc.ivw_fixed.heterogeneity", mr_ivw_fixed@Heter.Stat)
emit("mr.ldlc.ivw_random.estimate", mr_ivw_random@Estimate)
emit("mr.ldlc.ivw_random.standard_error", mr_ivw_random@StdError)
emit("mr.ldlc.egger.estimate", mr_egger_random@Estimate)
emit("mr.ldlc.egger.standard_error", mr_egger_random@StdError.Est)
emit("mr.ldlc.egger.intercept", mr_egger_random@Intercept)
emit("mr.ldlc.egger.intercept_standard_error", mr_egger_random@StdError.Int)
emit("mr.ldlc.egger.heterogeneity", mr_egger_random@Heter.Stat)
emit("mr.ldlc.weighted_median.estimate", mr_weighted_median@Estimate)
emit("mr.ldlc.weighted_median.standard_error", mr_weighted_median@StdError)

# Base-R inverse-variance meta-analysis and generalized DL meta-regression.
# These calculations are also directly reproducible with metafor::rma().
meta_y <- c(0.2, 0.5, 0.1, 0.7)
meta_se <- c(0.1, 0.2, 0.15, 0.25)
meta_v <- meta_se^2
meta_w <- 1 / meta_v
meta_mu <- sum(meta_w * meta_y) / sum(meta_w)
meta_q <- sum(meta_w * (meta_y - meta_mu)^2)
meta_c <- sum(meta_w) - sum(meta_w^2) / sum(meta_w)
meta_tau_dl <- max(0, (meta_q - 3) / meta_c)
meta_wr <- 1 / (meta_v + meta_tau_dl)
emit("meta.fixed.effect", meta_mu)
emit("meta.fixed.standard_error", sqrt(1 / sum(meta_w)))
emit("meta.fixed.cochran_q", meta_q)
emit("meta.fixed.cochran_p", pchisq(meta_q, 3, lower.tail = FALSE))
emit("meta.random_dl.tau_squared", meta_tau_dl)
emit("meta.random_dl.effect", sum(meta_wr * meta_y) / sum(meta_wr))
emit("meta.random_dl.standard_error", sqrt(1 / sum(meta_wr)))

meta_x <- cbind(1, c(-1, 0, 1, 2))
meta_b <- solve(t(meta_x) %*% (meta_w * meta_x),
                t(meta_x) %*% (meta_w * meta_y))
meta_residual <- meta_y - meta_x %*% meta_b
meta_qe <- sum(meta_w * meta_residual^2)
meta_information_inverse <- solve(t(meta_x) %*% (meta_w * meta_x))
meta_w2_information <- t(meta_x) %*% (meta_w^2 * meta_x)
meta_c_regression <- sum(meta_w) -
  sum(diag(meta_information_inverse %*% meta_w2_information))
meta_tau_regression <- max(0, (meta_qe - 2) / meta_c_regression)
meta_wr_regression <- 1 / (meta_v + meta_tau_regression)
meta_b_regression <- solve(t(meta_x) %*% (meta_wr_regression * meta_x),
                           t(meta_x) %*% (meta_wr_regression * meta_y))
meta_cov_regression <- solve(t(meta_x) %*% (meta_wr_regression * meta_x))
emit("meta.regression_dl.tau_squared", meta_tau_regression)
emit("meta.regression_dl.coefficients", meta_b_regression)
emit("meta.regression_dl.standard_errors", sqrt(diag(meta_cov_regression)))
