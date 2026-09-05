# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(TMB))

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else
  file.path("src", "test", "resources", "r-reference")
ids <- paste0("p", seq_len(60L))
sire <- c(rep(NA, 20L), rep(1:10, each = 4L))
dam <- c(rep(NA, 20L), rep(11:20, each = 4L))
q <- length(ids)
A <- matrix(0, q, q)
for (i in seq_len(q)) {
  A[i, i] <- if (is.na(sire[i]) || is.na(dam[i])) 1 else
    1 + 0.5 * A[sire[i], dam[i]]
  if (i > 1L) for (j in seq_len(i - 1L)) {
    A[i, j] <- 0.5 * (if (is.na(sire[i])) 0 else A[sire[i], j]) +
      0.5 * (if (is.na(dam[i])) 0 else A[dam[i], j])
    A[j, i] <- A[i, j]
  }
}
A_inverse <- as(solve(A), "dgCMatrix")
log_determinant_A_inverse <- as.numeric(determinant(A_inverse, logarithm = TRUE)$modulus)

set.seed(6102026)
observed_index <- rep(21:60, each = 30L)
n <- length(observed_index)
x <- rep(seq(-1, 1, length.out = 30L), 40L)
lower <- t(chol(A))
z1 <- rnorm(q)
z2 <- rnorm(q)
rho <- 0.55
a_count <- 0.7 * drop(lower %*% z1)
a_zero <- 1.0 * drop(lower %*% (rho * z1 + sqrt(1 - rho^2) * z2))
mu <- exp(0.4 + 0.3 * x + a_count[observed_index])
pi <- plogis(-1.0 - 0.2 * x + a_zero[observed_index])
y <- ifelse(runif(n) < pi, 0L, rpois(n, mu))
Z <- sparseMatrix(i = seq_len(n), j = observed_index, x = 1,
  dims = c(n, q))
X <- cbind(1, x)

temporary <- tempfile("jlinalg-tmb-")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
source_template <- file.path(root, "zero-inflated-pedigree-tmb.cpp")
template <- file.path(temporary, "zero_inflated_pedigree.cpp")
file.copy(source_template, template)
old <- setwd(temporary)
on.exit(setwd(old), add = TRUE)
TMB::compile(basename(template), flags = "-O2")
dyn.load(dynlib(sub("\\.cpp$", "", basename(template))))
objective <- MakeADFun(data = list(y = y, X_count = X, X_zero = X, Z = Z,
    A_inverse = A_inverse,
    log_determinant_A_inverse = log_determinant_A_inverse),
  parameters = list(beta_count = c(0, 0), beta_zero = c(-1, 0),
    a_count = rep(0, q), a_zero = rep(0, q),
    log_sd_count = log(0.5), log_sd_zero = log(0.7), fisher_z = 0),
  random = c("a_count", "a_zero"), DLL = "zero_inflated_pedigree",
  silent = TRUE)
fit <- nlminb(objective$par, objective$fn, objective$gr,
  control = list(iter.max = 1000, eval.max = 1500))
report <- sdreport(objective, par.fixed = fit$par)
setwd(old)
write.table(data.frame(y, x, id = ids[observed_index]),
  file.path(root, "zero-inflated-pedigree-tmb-data.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t")
values <- c(
  generator = "independent TMB 1.9.25 sparse-GMRF template under R 4.6.1",
  count_intercept = fit$par[[1L]],
  count_x = fit$par[[2L]],
  zero_intercept = fit$par[[3L]],
  zero_x = fit$par[[4L]],
  count_sd = exp(fit$par[["log_sd_count"]]),
  zero_sd = exp(fit$par[["log_sd_zero"]]),
  correlation = 0.99 * tanh(fit$par[["fisher_z"]]),
  log_likelihood = -fit$objective,
  convergence = fit$convergence)
connection <- file(file.path(root,
  "zero-inflated-pedigree-tmb.properties"), open = "wt")
on.exit(close(connection), add = TRUE)
for (name in names(values)) writeLines(paste0(name, "=", values[[name]]), connection)
