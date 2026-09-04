# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else
  file.path("src", "benchmark", "resources", "susie",
            "N3finemapping-reference.tsv")
local_library <- if (length(args) >= 2L) args[[2L]] else ".r-reference-lib"
if (dir.exists(local_library))
  .libPaths(c(normalizePath(local_library), .libPaths()))
if (!requireNamespace("susieR", quietly = TRUE))
  stop("install susieR before generating the reference")

data("N3finemapping", package = "susieR", envir = environment())
X <- N3finemapping$X
y <- N3finemapping$Y[, 1L]
prior_variance <- 0.2
fit <- susieR::susie(
  X, y, L = 10L, max_iter = 200L, tol = 1e-6,
  scaled_prior_variance = prior_variance / var(y),
  estimate_prior_variance = FALSE, estimate_residual_variance = TRUE,
  coverage = 0.95, min_abs_corr = 0.5, verbose = FALSE)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
connection <- file(output, open = "wt", encoding = "UTF-8")
writeLines(c(
  "# susieR 0.14.2; N3finemapping response replicate 1",
  sprintf("# residual_variance=%.17g", fit$sigma2),
  sprintf("# intercept=%.17g", fit$intercept),
  sprintf("# iterations=%d", fit$niter),
  sprintf("# objective=%.17g", tail(fit$elbo, 1L)),
  "index\tpip\tcoefficient"), connection)
coefficient <- coef(fit)[-1L]
for (index in seq_along(fit$pip)) {
  writeLines(sprintf("%d\t%.17g\t%.17g", index, fit$pip[[index]],
                     coefficient[[index]]), connection)
}
close(connection)
