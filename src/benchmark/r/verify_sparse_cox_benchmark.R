#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

arguments <- commandArgs(trailingOnly = TRUE)
options <- list(
  directory = "build/benchmarks/sparse-cox-validation",
  r_prefix = "r_validation",
  java_prefix = "java_validation",
  beta_tolerance = 0.000001,
  standard_error_tolerance = 0.0002,
  standardized_beta_tolerance = 0.00001,
  minimum_speedup = 1.0
)
if (length(arguments) %% 2L) stop("arguments must be --name value pairs")
for (index in seq(1L, length(arguments), by = 2L))
  options[[sub("^--", "", arguments[[index]])]] <- arguments[[index + 1L]]
for (name in c("beta_tolerance", "standard_error_tolerance",
    "standardized_beta_tolerance", "minimum_speedup"))
  options[[name]] <- as.numeric(options[[name]])

read_output <- function(prefix, suffix) fread(file.path(options$directory,
  paste0(prefix, suffix)), showProgress = FALSE)
r_results <- read_output(options$r_prefix, "_results.csv")
java_results <- read_output(options$java_prefix, "_results.csv")
r_timings <- read_output(options$r_prefix, "_timings.csv")
java_timings <- read_output(options$java_prefix, "_timings.csv")

keys <- c("model", "feature_key", "feature_id")
comparison <- merge(r_results, java_results, by = keys,
  suffixes = c("_r", "_java"))
if (!nrow(comparison) || nrow(comparison) != nrow(r_results) ||
    nrow(comparison) != nrow(java_results))
  stop("R and Java result keys do not match")
comparison[, `:=`(
  beta_absolute_difference = abs(beta_java - beta_r),
  standard_error_absolute_difference =
    abs(standard_error_java - standard_error_r),
  beta_difference_in_r_standard_errors =
    abs(beta_java - beta_r) / standard_error_r
)]

timing <- merge(
  r_timings[, .(r_median_seconds = median(seconds)), by = model],
  java_timings[, .(java_median_seconds = median(seconds)), by = model],
  by = "model")
timing[, speedup_r_over_java := r_median_seconds / java_median_seconds]
summary <- comparison[, .(
  comparisons = .N,
  max_beta_absolute_difference = max(beta_absolute_difference),
  max_standard_error_absolute_difference =
    max(standard_error_absolute_difference),
  max_beta_difference_in_r_standard_errors =
    max(beta_difference_in_r_standard_errors),
  all_java_converged = all(converged)), by = model]
summary <- merge(summary, timing, by = "model")
summary[, `:=`(
  beta_pass = max_beta_absolute_difference <= options$beta_tolerance,
  standard_error_pass = max_standard_error_absolute_difference <=
    options$standard_error_tolerance,
  standardized_beta_pass = max_beta_difference_in_r_standard_errors <=
    options$standardized_beta_tolerance,
  speed_pass = speedup_r_over_java > options$minimum_speedup
)]
summary[, passed := beta_pass & standard_error_pass &
  standardized_beta_pass & all_java_converged & speed_pass]

fwrite(comparison, file.path(options$directory,
  "sparse_cox_parity_details.csv"))
fwrite(summary, file.path(options$directory,
  "sparse_cox_validation_summary.csv"))
print(summary)
if (!all(summary$passed)) stop("sparse Cox R parity/performance gate failed")
