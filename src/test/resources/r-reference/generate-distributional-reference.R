# Deterministic independent references for distributional/vector additive fits.
local_library <- normalizePath(".r-reference-lib", mustWork = FALSE)
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(gamlss.dist))
suppressPackageStartupMessages(library(VGAM))

n <- 160L
i <- 0:(n - 1L)
x <- -1 + 2 * i / (n - 1)
noise <- (sin(1.7 * i) + 0.55 * cos(0.43 * i) + 0.25 * sin(2.31 * i)) /
  sqrt((1 + 0.55^2 + 0.25^2) / 2)
y <- 1.2 + 0.7 * x + exp(-0.4 + 0.3 * x) * noise
location_scale <- gamlss(y ~ x, sigma.formula = ~x, family = NO,
                         trace = FALSE, n.cyc = 200)
cat("gamlss.version=", as.character(packageVersion("gamlss")), "\n", sep = "")
cat("gamlss.mu=", paste(format(coef(location_scale, what = "mu"), digits = 17),
                         collapse = ","), "\n", sep = "")
cat("gamlss.sigma=", paste(format(coef(location_scale, what = "sigma"), digits = 17),
                            collapse = ","), "\n", sep = "")
cat("gamlss.logLikelihood=", format(as.numeric(logLik(location_scale)), digits = 17),
    "\n", sep = "")
cat("gamlss.mu.fitted=", paste(format(fitted(location_scale, what = "mu")[
    c(1L, 80L, 160L)], digits = 17), collapse = ","), "\n", sep = "")
cat("gamlss.sigma.fitted=", paste(format(fitted(location_scale, what = "sigma")[
    c(1L, 80L, 160L)], digits = 17), collapse = ","), "\n", sep = "")

n <- 180L
i <- 0:(n - 1L)
x <- -1 + 2 * i / (n - 1)
eta0 <- 0.2 + 0.7 * x
eta1 <- -0.3 - 0.5 * x
denominator <- 1 + exp(eta0) + exp(eta1)
p0 <- exp(eta0) / denominator
p1 <- exp(eta1) / denominator
state <- 7654321
uniform <- numeric(n)
for (row in seq_len(n)) {
  state <- (48271 * state) %% 2147483647
  uniform[row] <- state / 2147483647
}
y <- ifelse(uniform < p0, 0L, ifelse(uniform < p0 + p1, 1L, 2L))
multinomial <- vglm(factor(y, levels = 0:2) ~ x,
                    multinomial(refLevel = 3))
coefficient_matrix <- coef(multinomial, matrix = TRUE)
probability <- predict(multinomial, type = "response")
cat("VGAM.version=", as.character(packageVersion("VGAM")), "\n", sep = "")
cat("VGAM.coefficients=", paste(format(c(coefficient_matrix), digits = 17),
                                 collapse = ","), "\n", sep = "")
cat("VGAM.logLikelihood=", format(as.numeric(logLik(multinomial)), digits = 17),
    "\n", sep = "")
cat("VGAM.probability=", paste(format(probability[cbind(
    c(1L, 90L, 180L), c(1L, 2L, 3L))], digits = 17), collapse = ","),
    "\n", sep = "")
