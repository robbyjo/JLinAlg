#!/usr/bin/env Rscript

options(digits = 17)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1")

options <- list(exposures = 300L, outcomes = 150L, instruments = 10L,
  warmups = 2L, measurements = 5L,
  output = file.path("build", "benchmarks", "xwas-mr"))
args <- commandArgs(trailingOnly = TRUE)
index <- 1L
while (index <= length(args)) {
  if (index == length(args)) stop("missing value for ", args[[index]])
  key <- args[[index]]
  value <- args[[index + 1L]]
  if (key == "--exposures") options$exposures <- as.integer(value)
  else if (key == "--outcomes") options$outcomes <- as.integer(value)
  else if (key == "--instruments") options$instruments <- as.integer(value)
  else if (key == "--warmups") options$warmups <- as.integer(value)
  else if (key == "--measurements") options$measurements <- as.integer(value)
  else if (key == "--output") options$output <- value
  else stop("unknown option: ", key)
  index <- index + 2L
}
if (options$exposures < 1L || options$outcomes < 1L ||
    options$instruments < 3L || options$warmups < 0L ||
    options$measurements < 1L) stop("invalid benchmark dimensions")
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)

exposure_effect <- function(exposure, instrument) {
  0.05 + 0.004 * instrument + 0.0001 * (exposure %% 17)
}
outcome_effect <- function(exposure, outcome, instrument) {
  causal <- ((exposure * 13 + outcome * 7) %% 11 - 5) * 0.02
  noise <- 0.01 * sin((exposure + 1) * 0.17 + (outcome + 1) * 0.11 +
    (instrument + 1) * 0.37)
  causal * exposure_effect(exposure, instrument) + noise
}

instrument_index <- 0:(options$instruments - 1L)
x <- matrix(NA_real_, options$exposures, options$instruments)
y <- array(NA_real_, c(options$exposures, options$outcomes,
  options$instruments))
for (exposure in 0:(options$exposures - 1L)) {
  x[exposure + 1L, ] <- exposure_effect(exposure, instrument_index)
  for (outcome in 0:(options$outcomes - 1L))
    y[exposure + 1L, outcome + 1L, ] <- outcome_effect(
      exposure, outcome, instrument_index)
}
se_y <- 0.04 + 0.001 * (instrument_index %% 4)
weight <- 1 / se_y^2
degrees_of_freedom <- options$instruments - 1L

fit_one <- function(x_value, y_value) {
  information <- sum(weight * x_value^2)
  beta <- sum(weight * x_value * y_value) / information
  q <- sum(((y_value - beta * x_value) / se_y)^2)
  dispersion <- max(1, q / degrees_of_freedom)
  se <- sqrt(dispersion / information)
  statistic <- beta / se
  p <- min(1, 2 * pnorm(abs(statistic), lower.tail = FALSE))
  c(beta = beta, se = se, statistic = statistic, p_value = p,
    cochran_q = q,
    heterogeneity_p_value = pchisq(q, degrees_of_freedom,
      lower.tail = FALSE),
    i_squared = if (q > 0) max(0, (q - degrees_of_freedom) / q) else 0)
}

scan <- function() {
  checksum <- 0
  for (exposure in seq_len(options$exposures)) {
    x_value <- x[exposure, ]
    for (outcome in seq_len(options$outcomes)) {
      value <- fit_one(x_value, y[exposure, outcome, ])
      checksum <- checksum + value[["beta"]] + value[["p_value"]]
    }
  }
  checksum
}

if (options$warmups > 0L)
  for (warmup in seq_len(options$warmups)) invisible(scan())
timings <- vector("list", options$measurements)
for (measurement in seq_len(options$measurements)) {
  gc()
  started <- proc.time()[["elapsed"]]
  checksum <- scan()
  seconds <- proc.time()[["elapsed"]] - started
  if (!is.finite(checksum)) stop("non-finite benchmark checksum")
  pairs <- options$exposures * options$outcomes
  timings[[measurement]] <- data.frame(runtime = "R-base-loop",
    measurement = measurement, exposures = options$exposures,
    outcomes = options$outcomes, pairs = pairs,
    instruments = options$instruments, threads = 1L, seconds = seconds,
    pairs_per_second = pairs / seconds)
  cat(sprintf(paste0("R-base-loop measurement=%d exposures=%d outcomes=%d ",
    "pairs=%d instruments=%d threads=1 seconds=%.9f ",
    "pairs_per_second=%.3f\n"), measurement, options$exposures,
    options$outcomes, pairs, options$instruments, seconds, pairs / seconds))
}
write.csv(do.call(rbind, timings),
  file.path(options$output, "r_timings.csv"), row.names = FALSE, quote = FALSE)

validation <- list()
row <- 0L
for (exposure in 0:(min(5L, options$exposures) - 1L))
  for (outcome in 0:(min(5L, options$outcomes) - 1L)) {
    row <- row + 1L
    value <- fit_one(x[exposure + 1L, ], y[exposure + 1L, outcome + 1L, ])
    validation[[row]] <- data.frame(exposure = exposure, outcome = outcome,
      beta = value[["beta"]], se = value[["se"]],
      statistic = value[["statistic"]], p_value = value[["p_value"]],
      cochran_q = value[["cochran_q"]],
      heterogeneity_p_value = value[["heterogeneity_p_value"]],
      i_squared = value[["i_squared"]])
  }
write.csv(do.call(rbind, validation),
  file.path(options$output, "r_validation.csv"), row.names = FALSE,
  quote = FALSE)
