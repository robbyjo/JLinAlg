# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

local_reference_library <- file.path(getwd(), ".r-reference-lib")
if (dir.exists(local_reference_library)) {
  .libPaths(c(local_reference_library, .libPaths()))
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glmnet))

parse_arguments <- function(arguments) {
  result <- list(prepared_dir = "build/benchmarks/topmed100", genes = 100L,
    measurements = 3L, warmups = 5L, repetitions = 10L, lambdas = 100L,
    minimum_ratio = 1e-4,
    tolerance = 1e-8, models = "ridge,lasso,elastic-net",
    output_prefix = "build/benchmarks/topmed100/r_penalized")
  index <- 1L
  while (index <= length(arguments)) {
    key <- gsub("-", "_", sub("^--", "", arguments[[index]]), fixed = TRUE)
    if (index == length(arguments)) stop("missing value for --", key)
    result[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  result$genes <- as.integer(result$genes)
  result$measurements <- as.integer(result$measurements)
  result$warmups <- as.integer(result$warmups)
  result$repetitions <- as.integer(result$repetitions)
  result$lambdas <- as.integer(result$lambdas)
  result$minimum_ratio <- as.numeric(result$minimum_ratio)
  result$tolerance <- as.numeric(result$tolerance)
  result
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
analysis <- fread(file.path(options$prepared_dir, "analysis.csv"),
  showProgress = FALSE)
features <- head(fread(file.path(options$prepared_dir, "features.csv"),
  colClasses = "character", showProgress = FALSE), options$genes)
covariates <- c("Sex", "Age", "WBC_Pred", "LY_PER_Pred", "MO_PER_Pred",
  "EO_PER_Pred", "BA_PER_Pred")
terms <- c(features$feature_key, covariates)
x <- as.matrix(analysis[, ..terms])
storage.mode(x) <- "double"
y <- analysis$BMI
x_means <- colMeans(x)
x_scales <- sqrt(colMeans(sweep(x, 2L, x_means)^2))
x_fit <- sweep(sweep(x, 2L, x_means), 2L, x_scales, "/")
y_mean <- mean(y)
y_scale <- sqrt(mean((y - y_mean)^2))
y_fit <- (y - y_mean) / y_scale
penalty_factor <- c(rep(1, nrow(features)), rep(0, length(covariates)))
model_names <- strsplit(options$models, ",", fixed = TRUE)[[1L]]
alphas <- c(ridge = 0, lasso = 1, `elastic-net` = 0.5)

fit_model <- function(alpha) glmnet(x_fit, y_fit, family = "gaussian", alpha = alpha,
  nlambda = options$lambdas, lambda.min.ratio = options$minimum_ratio,
  penalty.factor = penalty_factor, standardize = FALSE, intercept = FALSE,
  control = list(thresh = options$tolerance, fdev = 0, devmax = 1))

timings <- list()
coefficients <- list()
for (model in model_names) {
  alpha <- alphas[[model]]
  if (is.null(alpha)) stop("unknown model: ", model)
  for (warmup in seq_len(options$warmups)) invisible(fit_model(alpha))
  for (measurement in seq_len(options$measurements)) {
    invisible(gc())
    elapsed <- system.time(for (repetition in seq_len(options$repetitions)) {
      fit <- fit_model(alpha)
    })[["elapsed"]] / options$repetitions
    timings[[length(timings) + 1L]] <- data.table(runtime = "R", model,
      threads = 1L, measurement, samples = nrow(x), genes = nrow(features),
      lambdas = length(fit$lambda), seconds = elapsed,
      lambda_fits_per_second = length(fit$lambda) / elapsed)
    if (measurement == 1L) {
      beta <- as.matrix(coef(fit, s = tail(fit$lambda, 1L)))
      beta[-1L, 1L] <- beta[-1L, 1L] * y_scale / x_scales
      beta[1L, 1L] <- y_mean - sum(x_means * beta[-1L, 1L])
      coefficients[[length(coefficients) + 1L]] <- data.table(runtime = "R",
        model, alpha, lambda = tail(fit$lambda, 1L) * y_scale, term = rownames(beta),
        beta = beta[, 1L])
    }
    cat(sprintf("R model=%s threads=1 measurement=%d seconds=%.6f lambda_fits_per_second=%.1f\n",
      model, measurement, elapsed, length(fit$lambda) / elapsed))
  }
}
for (warmup in seq_len(options$warmups)) {
  invisible(lapply(model_names, function(model) fit_model(alphas[[model]])))
}
for (measurement in seq_len(options$measurements)) {
  invisible(gc())
  elapsed <- system.time(for (repetition in seq_len(options$repetitions)) {
    fits <- lapply(model_names, function(model) fit_model(alphas[[model]]))
  })[["elapsed"]] / options$repetitions
  lambda_count <- sum(vapply(fits, function(fit) length(fit$lambda), 0L))
  timings[[length(timings) + 1L]] <- data.table(runtime = "R", model = "all",
    threads = 1L, measurement, samples = nrow(x), genes = nrow(features),
    lambdas = lambda_count, seconds = elapsed,
    lambda_fits_per_second = lambda_count / elapsed)
  cat(sprintf("R model=all threads=1 measurement=%d seconds=%.6f lambda_fits_per_second=%.1f\n",
    measurement, elapsed, lambda_count / elapsed))
}
fwrite(rbindlist(timings), paste0(options$output_prefix, "_timings.csv"))
fwrite(rbindlist(coefficients), paste0(options$output_prefix, "_coefficients.csv"))
