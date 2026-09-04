# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

args <- commandArgs(trailingOnly = TRUE)
warmups <- if (length(args) >= 1L) as.integer(args[[1L]]) else 3L
measurements <- if (length(args) >= 2L) as.integer(args[[2L]]) else 7L
local_library <- if (length(args) >= 3L) args[[3L]] else ".r-reference-lib"
if (dir.exists(local_library))
  .libPaths(c(normalizePath(local_library), .libPaths()))
if (!requireNamespace("susieR", quietly = TRUE))
  stop("install susieR before running this benchmark")
if (is.na(warmups) || warmups < 1L ||
    is.na(measurements) || measurements < 1L)
  stop("warmups and measurements must be positive integers")

data("N3finemapping", package = "susieR", envir = environment())
X <- N3finemapping$X
y <- N3finemapping$Y[, 1L]
fit <- function() suppressMessages(susieR::susie(
  X, y, L = 10L, max_iter = 200L, tol = 1e-6,
  scaled_prior_variance = 0.2 / var(y), estimate_prior_variance = FALSE,
  estimate_residual_variance = TRUE, coverage = 0.95,
  min_abs_corr = 0.5, verbose = FALSE))

for (iteration in seq_len(warmups)) result <- fit()
seconds <- numeric(measurements)
for (iteration in seq_len(measurements)) {
  started <- proc.time()[["elapsed"]]
  result <- fit()
  seconds[[iteration]] <- proc.time()[["elapsed"]] - started
}
checksum <- result$sigma2 + result$intercept +
  sum(result$pip * seq_along(result$pip)) + sum(coef(result)[-1L])
cat("runtime,dataset,observations,variables,effects,warmups,measurements,median_seconds,iterations,converged,residual_variance,credible_sets,backend,checksum\n")
cat(paste("R", "N3finemapping", nrow(X), ncol(X), nrow(result$alpha),
  warmups, measurements, sprintf("%.9f", median(seconds)), result$niter,
  result$converged, sprintf("%.12g", result$sigma2), length(result$sets$cs),
  encodeString("R BLAS", quote = '"'), sprintf("%.12g", checksum),
  sep = ","), "\n", sep = "")
