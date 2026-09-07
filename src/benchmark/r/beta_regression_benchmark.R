#!/usr/bin/env Rscript

benchmark_library <- Sys.getenv("JLINALG_BETA_R_LIBRARY", "")
if (nzchar(benchmark_library)) .libPaths(c(benchmark_library, .libPaths()))
suppressPackageStartupMessages(library(betareg))

rows <- as.integer(Sys.getenv("JLINALG_BETA_ROWS", "100000"))
repetitions <- as.integer(Sys.getenv("JLINALG_BETA_REPETITIONS", "5"))
index <- 0:(rows - 1L)
x1 <- -1 + 2 * index / max(1, rows - 1)
x2 <- sin(0.017 * index)
x3 <- cos(0.031 * index)
mu <- plogis(-0.4 + 0.8 * x1 - 0.35 * x2 + 0.2 * x3)
phi <- exp(3.8 + 0.25 * x1)
standardized <- (sin(1.73 * index + 0.2) +
  0.55 * cos(0.47 * index)) / 1.14
response <- pmax(1e-5, pmin(1 - 1e-5,
  mu + sqrt(mu * (1 - mu) / (phi + 1)) * standardized))
data <- data.frame(response, x1, x2, x3)

fit_once <- function() {
  betareg(response ~ x1 + x2 + x3 | x1, data = data)
}

invisible(fit_once())
times <- numeric(repetitions)
for (repetition in seq_len(repetitions)) {
  gc()
  timing <- system.time(fit <- fit_once())
  times[[repetition]] <- unname(timing[["elapsed"]]) * 1000
}

cat(sprintf("betareg_version=%s rows=%d repetitions=%d\n",
  packageVersion("betareg"), rows, repetitions))
cat(sprintf("median_ms=%.3f\n", median(times)))
cat("coefficients=", paste(format(coef(fit), digits = 16), collapse = ","),
  "\n", sep = "")
