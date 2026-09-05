# Deterministic conventional survival::coxph.fit scan comparator.
suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly = TRUE)
rows <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
variables <- if (length(args) >= 2) as.integer(args[[2]]) else 512L
measurements <- if (length(args) >= 3) as.integer(args[[3]]) else 3L
set.seed(20260905)
covariates <- cbind(rnorm(rows), rnorm(rows))
predictors <- matrix(rnorm(rows * variables), rows, variables)
time <- 0.25 + 10 * runif(rows)
event <- runif(rows) < 0.55
event[[1]] <- TRUE
response <- Surv(time, event)
strata <- rep.int(0L, rows)
row_names <- as.character(seq_len(rows))

scan <- function() {
  last <- 0
  for (j in seq_len(variables)) {
    fit <- coxph.fit(cbind(covariates, predictors[, j]), response,
      strata = strata, method = "efron", rownames = row_names,
      control = coxph.control())
    last <- fit$coefficients[[3]]
  }
  invisible(last)
}

scan() # warm-up
elapsed <- replicate(measurements, system.time(scan())[["elapsed"]])
median_seconds <- median(elapsed)
cat("runtime,rows,variables,median_seconds,variables_per_second\n")
cat(sprintf("R_survival_full_refit,%d,%d,%.6f,%.2f\n", rows, variables,
  median_seconds, variables / median_seconds))
