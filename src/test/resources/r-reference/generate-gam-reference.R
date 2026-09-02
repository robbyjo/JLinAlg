# Deterministic mgcv reference for JLinAlg's GAM validation suite.
#
# This script is intentionally separate from generate-reference.R so the core
# reference suite does not require mgcv until GAM validation is requested.
# Run from the repository root with:
#   Rscript src/test/resources/r-reference/generate-gam-reference.R

options(digits = 17)
suppressPackageStartupMessages(library(mgcv))

emit <- function(name, value) {
  cat(name, "=")
  dput(unname(value))
}

observations <- 80L
x <- seq(0, 1, length.out = observations)
row <- seq_len(observations) - 1L
y <- 1.5 + sin(2 * pi * x) + 0.08 * sin(17 * row)

# P-spline and derivative order mirror PSplineTerm.of(name, x, 10).
fit <- gam(y ~ s(x, bs = "ps", k = 10, m = 2), method = "REML")
summary_fit <- summary(fit)

cat("R.version=", R.version.string, "\n", sep = "")
cat("mgcv.version=", as.character(packageVersion("mgcv")), "\n", sep = "")
emit("gam.ps.x", x)
emit("gam.ps.y", y)
emit("gam.ps.coefficients", coef(fit))
emit("gam.ps.smoothing_parameter", fit$sp)
emit("gam.ps.edf", summary_fit$s.table[, "edf"])
emit("gam.ps.reference_df", summary_fit$s.table[, "Ref.df"])
emit("gam.ps.statistic", summary_fit$s.table[, "F"])
emit("gam.ps.p_value", summary_fit$s.table[, "p-value"])
emit("gam.ps.fitted", fitted(fit))
emit("gam.ps.deviance", deviance(fit))
emit("gam.ps.log_likelihood", as.numeric(logLik(fit)))
