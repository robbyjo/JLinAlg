# Optional cross-package compatibility checks for JLinAlg additive stages.
# Run with: Rscript generate-additive-optional-references.R

local_library <- normalizePath(".r-reference-lib", mustWork = FALSE)
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

set.seed(20260902)
n <- 120L
x <- seq(0, 1, length.out = n)
group <- factor(rep(seq_len(12L), each = 10L))
y <- 1.2 + sin(2 * pi * x) + rnorm(n, sd = exp(-1 + 0.4 * x))
count <- rpois(n, exp(0.3 + 0.5 * sin(2 * pi * x)))

run_if_available <- function(package, expression) {
  if (!requireNamespace(package, quietly = TRUE)) {
    cat("SKIP ", package, ": package is not installed\n", sep = "")
    return(invisible(NULL))
  }
  cat("BEGIN ", package, " ", as.character(packageVersion(package)), "\n", sep = "")
  force(expression)
  cat("END ", package, "\n", sep = "")
}

run_if_available("mgcv", {
  fit <- mgcv::gam(y ~ s(x, bs = "ps", k = 10), method = "REML")
  cat("edf=", summary(fit)$s.table[1, "edf"],
      " logLik=", as.numeric(logLik(fit)), "\n", sep = "")
})

run_if_available("gamm4", {
  fit <- gamm4::gamm4(count ~ s(x, k = 8), random = ~(1 | group),
                      family = poisson())
  cat("beta=", lme4::fixef(fit$mer)[1], "\n", sep = "")
})

run_if_available("gamlss", {
  fit <- gamlss::gamlss(y ~ gamlss::pb(x), sigma.formula = ~x,
                        family = gamlss.dist::NO, trace = FALSE)
  cat("mu.df=", fit$mu.df, " sigma.df=", fit$sigma.df,
      " globalDeviance=", fit$G.deviance, "\n", sep = "")
})

run_if_available("gam", {
  fit <- gam::gam(y ~ gam::s(x, df = 6))
  cat("deviance=", deviance(fit), "\n", sep = "")
})

run_if_available("VGAM", {
  category <- ordered(cut(y, breaks = quantile(y, 0:4 / 4),
                          include.lowest = TRUE))
  fit <- VGAM::vglm(category ~ x, VGAM::acat(parallel = TRUE))
  cat("logLik=", as.numeric(logLik(fit)), "\n", sep = "")
})
