#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1")
setDTthreads(1L)

options <- list(
  prepared_dir = "build/benchmarks/topmed-meta-analysis",
  measurements = 3L,
  models = c("fixed", "reml"),
  maximum_iterations = 200L,
  tolerance = 1e-10)
args <- commandArgs(trailingOnly = TRUE)
index <- 1L
while (index <= length(args)) {
  key <- args[[index]]
  if (index == length(args)) stop("missing value for ", key)
  value <- args[[index + 1L]]
  if (key == "--prepared-dir") options$prepared_dir <- value
  else if (key == "--measurements") options$measurements <- as.integer(value)
  else if (key == "--models")
    options$models <- strsplit(value, ",", fixed = TRUE)[[1L]]
  else stop("unknown option: ", key)
  index <- index + 2L
}

analysis <- fread(file.path(options$prepared_dir, "analysis.csv"))
effect_columns <- grep("^effect_", names(analysis), value = TRUE)
error_columns <- grep("^se_", names(analysis), value = TRUE)
effects <- as.matrix(analysis[, ..effect_columns])
variances <- as.matrix(analysis[, ..error_columns])^2
studies <- ncol(effects)
degrees_of_freedom <- studies - 1
critical <- qnorm(0.975)

weighted_sums <- function(tau_squared) {
  total_variance <- variances + tau_squared
  weights <- 1 / total_variance
  sum_weight <- rowSums(weights)
  mean <- rowSums(weights * effects) / sum_weight
  residual <- effects - mean
  q <- rowSums(weights * residual^2)
  list(sum_weight = sum_weight, mean = mean, q = q,
    total_variance = total_variance)
}

restricted_objective <- function(tau_squared) {
  value <- weighted_sums(tau_squared)
  rowSums(log(value$total_variance)) + log(value$sum_weight) + value$q +
    degrees_of_freedom * log(2 * pi)
}

finish_fit <- function(tau_squared, fixed, random_effects) {
  fitted <- if (all(tau_squared == 0)) fixed else weighted_sums(tau_squared)
  standard_error <- sqrt(1 / fitted$sum_weight)
  statistic <- fitted$mean / standard_error
  p_value <- pmin(1, 2 * pnorm(abs(statistic), lower.tail = FALSE))
  prediction_error <- if (random_effects)
    sqrt(tau_squared + standard_error^2) else rep(NA_real_, length(tau_squared))
  list(beta = fitted$mean, se = standard_error, statistic = statistic,
    p_value = p_value,
    confidence_lower = fitted$mean - critical * standard_error,
    confidence_upper = fitted$mean + critical * standard_error,
    prediction_lower = fitted$mean - critical * prediction_error,
    prediction_upper = fitted$mean + critical * prediction_error,
    cochran_q = fixed$q,
    cochran_q_p_value = pchisq(fixed$q, degrees_of_freedom,
      lower.tail = FALSE),
    tau_squared = tau_squared,
    i_squared = ifelse(fixed$q > 0,
      100 * pmax(0, (fixed$q - degrees_of_freedom) / fixed$q), 0),
    h_squared = pmax(1, fixed$q / degrees_of_freedom))
}

fit_fixed <- function() {
  fixed <- weighted_sums(rep(0, nrow(effects)))
  finish_fit(rep(0, nrow(effects)), fixed, FALSE)
}

fit_reml <- function() {
  fixed <- weighted_sums(rep(0, nrow(effects)))
  means <- rowMeans(effects)
  upper <- pmax(1e-8,
    rowSums((effects - means)^2) / degrees_of_freedom)
  at_zero <- restricted_objective(rep(0, nrow(effects)))
  previous <- at_zero
  at_upper <- restricted_objective(upper)
  repeat {
    grow <- at_upper < previous & upper < 1e12
    if (!any(grow)) break
    previous[grow] <- at_upper[grow]
    upper[grow] <- upper[grow] * 4
    candidate <- restricted_objective(upper)
    at_upper[grow] <- candidate[grow]
  }
  left <- rep(0, length(upper))
  right <- upper
  ratio <- (sqrt(5) - 1) / 2
  first <- right - ratio * (right - left)
  second <- left + ratio * (right - left)
  first_value <- restricted_objective(first)
  second_value <- restricted_objective(second)
  for (iteration in seq_len(options$maximum_iterations)) {
    pending <- right - left > options$tolerance * pmax(1, right)
    if (!any(pending)) break
    choose_first <- pending & first_value < second_value
    choose_second <- pending & !choose_first
    right[choose_first] <- second[choose_first]
    second[choose_first] <- first[choose_first]
    second_value[choose_first] <- first_value[choose_first]
    first[choose_first] <- right[choose_first] -
      ratio * (right[choose_first] - left[choose_first])
    left[choose_second] <- first[choose_second]
    first[choose_second] <- second[choose_second]
    first_value[choose_second] <- second_value[choose_second]
    second[choose_second] <- left[choose_second] +
      ratio * (right[choose_second] - left[choose_second])
    next_point <- ifelse(choose_first, first, second)
    next_value <- restricted_objective(next_point)
    first_value[choose_first] <- next_value[choose_first]
    second_value[choose_second] <- next_value[choose_second]
  }
  tau_squared <- 0.5 * (left + right)
  at_candidate <- restricted_objective(tau_squared)
  tau_squared[at_zero <= at_candidate] <- 0
  finish_fit(tau_squared, fixed, TRUE)
}

fit_model <- function(model) {
  if (model == "fixed") fit_fixed()
  else if (model == "reml") fit_reml()
  else stop("unknown model: ", model)
}

timings <- list()
results <- list()
timing_index <- 0L
for (model in options$models) {
  invisible(fit_model(model))
  first <- NULL
  for (measurement in seq_len(options$measurements)) {
    gc()
    started <- proc.time()[["elapsed"]]
    value <- fit_model(model)
    seconds <- proc.time()[["elapsed"]] - started
    if (is.null(first)) first <- value
    if (!all(is.finite(value$p_value)) ||
        !all(is.finite(value$tau_squared)))
      stop("non-finite benchmark result")
    timing_index <- timing_index + 1L
    timings[[timing_index]] <- data.table(runtime = "R-base-vectorized",
      engine = "vectorized", model = model,
      tau_estimator = if (model == "fixed") "REML" else "REML",
      threads = 1L, measurement = measurement,
      analyses = nrow(effects), seconds = seconds,
      analyses_per_second = nrow(effects) / seconds)
    cat(sprintf(paste0("R-base-vectorized model=%s threads=1 ",
      "measurement=%d analyses=%d seconds=%.6f ",
      "analyses_per_second=%.3f\n"), model, measurement,
      nrow(effects), seconds, nrow(effects) / seconds))
  }
  results[[model]] <- data.table(runtime = "R-base-vectorized",
    model = model, feature = analysis$feature, beta = first$beta,
    se = first$se, p_value = first$p_value,
    cochran_q = first$cochran_q,
    cochran_q_p_value = first$cochran_q_p_value,
    tau_squared = first$tau_squared, i_squared = first$i_squared)
}

fwrite(rbindlist(timings),
  file.path(options$prepared_dir, "r_timings.csv"))
fwrite(rbindlist(results),
  file.path(options$prepared_dir, "r_results.csv"))
