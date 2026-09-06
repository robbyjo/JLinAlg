# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later
# Export identical model matrices and benchmark glmmTMB's Salamanders examples.
suppressPackageStartupMessages(library(glmmTMB))
data(Salamanders, package = "glmmTMB")
options(glmmTMB.cores = 1L)

zip_call <- function() glmmTMB(
  count ~ mined + (1 | site), ziformula = ~ mined,
  family = poisson, data = Salamanders)
zinb_call <- function() glmmTMB(
  count ~ spp + mined + (1 | site), ziformula = ~ spp + mined,
  family = nbinom2, data = Salamanders)
zip_call_nose <- function() glmmTMB(
  count ~ mined + (1 | site), ziformula = ~ mined,
  family = poisson, data = Salamanders, se = FALSE)
zinb_call_nose <- function() glmmTMB(
  count ~ spp + mined + (1 | site), ziformula = ~ spp + mined,
  family = nbinom2, data = Salamanders, se = FALSE)

zip <- zip_call()
zinb <- zinb_call()

timings <- function(fun, warmup = 2L, runs = 7L) {
  invisible(replicate(warmup, fun(), simplify = FALSE))
  elapsed <- replicate(runs, system.time(fun())[["elapsed"]])
  c(median = median(elapsed), min = min(elapsed), max = max(elapsed))
}

zip_time <- timings(zip_call)
zinb_time <- timings(zinb_call)
zip_time_nose <- timings(zip_call_nose)
zinb_time_nose <- timings(zinb_call_nose)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else file.path("build", "tmp", "salamander-comparison")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(
    count = Salamanders$count,
    site = as.character(Salamanders$site),
    model.matrix(~ mined, Salamanders),
    check.names = FALSE),
  file.path(root, "zip.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(
    count = Salamanders$count,
    site = as.character(Salamanders$site),
    model.matrix(~ spp + mined, Salamanders),
    check.names = FALSE),
  file.path(root, "zinb.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

emit <- function(label, fit, timing) {
  cat(label, "_count=", paste(format(fixef(fit)$cond, digits = 17), collapse = ","), "\n", sep = "")
  cat(label, "_zero=", paste(format(fixef(fit)$zi, digits = 17), collapse = ","), "\n", sep = "")
  cat(label, "_count_sd=", format(attr(VarCorr(fit)$cond$site, "stddev")[[1L]], digits = 17), "\n", sep = "")
  if (label == "zinb") cat(label, "_size=", format(sigma(fit), digits = 17), "\n", sep = "")
  cat(label, "_loglik=", format(as.numeric(logLik(fit)), digits = 17), "\n", sep = "")
  cat(label, "_se_count=", paste(format(summary(fit)$coefficients$cond[, 2L], digits = 17), collapse = ","), "\n", sep = "")
  cat(label, "_se_zero=", paste(format(summary(fit)$coefficients$zi[, 2L], digits = 17), collapse = ","), "\n", sep = "")
  cat(label, "_time_median=", timing[["median"]], "\n", sep = "")
  cat(label, "_time_min=", timing[["min"]], "\n", sep = "")
  cat(label, "_time_max=", timing[["max"]], "\n", sep = "")
}
cat("r_threads=1\n")
cat("r_version=", R.version.string, "\n", sep = "")
cat("glmmtmb_version=", as.character(packageVersion("glmmTMB")), "\n", sep = "")
cat("rows=", nrow(Salamanders), "\n", sep = "")
emit("zip", zip, zip_time)
emit("zinb", zinb, zinb_time)
cat("zip_time_no_se_median=", zip_time_nose[["median"]], "\n", sep = "")
cat("zinb_time_no_se_median=", zinb_time_nose[["median"]], "\n", sep = "")
