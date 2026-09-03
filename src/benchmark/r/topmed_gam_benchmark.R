# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mgcv))
setDTthreads(1L)

options <- list(
  prepared_dir = "build/benchmarks/topmed100",
  genes = 100L,
  measurements = 1L
)
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) %% 2L) stop("arguments must be --name value pairs")
for (index in seq(1L, length(arguments), by = 2L)) {
  key <- gsub("-", "_", sub("^--", "", arguments[[index]]), fixed = TRUE)
  options[[key]] <- arguments[[index + 1L]]
}
options$genes <- as.integer(options$genes)
options$measurements <- as.integer(options$measurements)

analysis <- fread(
  file.path(options$prepared_dir, "analysis.csv"),
  showProgress = FALSE
)
features <- head(fread(
  file.path(options$prepared_dir, "features.csv"),
  colClasses = "character", showProgress = FALSE
), options$genes)
fixed_rhs <- paste(
  "Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred +",
  "EO_PER_Pred + BA_PER_Pred"
)
fit_gene <- function(feature) {
  formula <- as.formula(paste(
    "BMI ~ s(", feature, ", bs='ps', k=10, m=2) +", fixed_rhs
  ))
  gam(formula, data = analysis, method = "REML")
}

invisible(fit_gene(features$feature_key[[1L]]))
timing_rows <- list()
result_rows <- NULL
for (measurement in seq_len(options$measurements)) {
  invisible(gc())
  elapsed <- system.time({
    fits <- lapply(features$feature_key, fit_gene)
  })[["elapsed"]]
  timing_rows[[measurement]] <- data.table(
    runtime = "R", model = "gam", backend = "mgcv",
    threads = 1L, measurement, genes = nrow(features),
    seconds = elapsed, genes_per_second = nrow(features) / elapsed
  )
  if (measurement == 1L) {
    result_rows <- rbindlist(lapply(seq_along(fits), function(index) {
      fit <- fits[[index]]
      data.table(
        runtime = "R", model = "gam",
        feature_key = features$feature_key[[index]],
        feature_id = features$feature_id[[index]],
        edf = unname(summary(fit)$s.table[1L, "edf"]),
        smoothing_parameter = unname(fit$sp[[1L]]),
        log_likelihood = as.numeric(logLik(fit)),
        fitted_checksum = sum(fitted(fit)),
        residual_sum_squares = sum(residuals(fit)^2)
      )
    }))
  }
}
timings <- rbindlist(timing_rows)
fwrite(timings, file.path(options$prepared_dir, "r_gam_timings.csv"))
fwrite(result_rows, file.path(options$prepared_dir, "r_gam_results.csv"))
print(timings)
