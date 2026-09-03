# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else
  file.path("src", "benchmark", "resources", "timeseries")
warmups <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1L
measurements <- if (length(args) >= 3L) as.integer(args[[3L]]) else 3L
repetitions <- if (length(args) >= 4L) as.integer(args[[4L]]) else 50L
selected <- if (length(args) >= 5L && nzchar(args[[5L]]))
  strsplit(args[[5L]], ",", fixed = TRUE)[[1L]] else character()

if (is.na(warmups) || warmups < 1L || is.na(measurements) || measurements < 1L ||
    is.na(repetitions) || repetitions < 1L)
  stop("warmups, measurements, and repetitions must be positive integers")

read_series <- function(name, expected) {
  path <- file.path(data_dir, paste0(name, ".csv"))
  frame <- read.csv(path, check.names = FALSE)
  if (!identical(names(frame), c("rownames", "time", "value")) ||
      nrow(frame) != expected || any(!is.finite(frame$value)) ||
      any(diff(frame$time) <= 0))
    stop("invalid benchmark data: ", path)
  frame$value
}

data <- list(
  AirPassengers = read_series("AirPassengers", 144L),
  Nile = read_series("Nile", 100L),
  nottem = read_series("nottem", 240L),
  sunspots = read_series("sunspots", 2820L),
  UKgas = read_series("UKgas", 108L),
  WWWusage = read_series("WWWusage", 100L)
)

fit <- function(values, order, seasonal = NULL) {
  arguments <- list(x = values, order = order, method = "CSS",
                    optim.control = list(maxit = 5000L, reltol = 1e-8))
  if (!is.null(seasonal)) arguments$seasonal <- seasonal
  do.call(stats::arima, arguments)
}

fit_ml <- function(values, order) {
  stats::arima(values, order = order, method = "ML",
               optim.control = list(maxit = 5000L, reltol = 1e-8))
}

cases <- list(
  conditional_ar2_nile = list(
    dataset = "Nile", model = "AR(2)", period = 1L,
    run = function() fit(data$Nile, c(2L, 0L, 0L))),
  conditional_arima310_wwwusage = list(
    dataset = "WWWusage", model = "ARIMA(3,1,0)", period = 1L,
    run = function() fit(data$WWWusage, c(3L, 1L, 0L))),
  conditional_sarima_airline = list(
    dataset = "AirPassengers", model = "log SARIMA(0,1,1)(0,1,1)[12]",
    period = 12L,
    run = function() fit(log(data$AirPassengers), c(0L, 1L, 1L),
      list(order = c(0L, 1L, 1L), period = 12L))),
  conditional_sarima_ukgas = list(
    dataset = "UKgas", model = "log SARIMA(1,1,1)(1,1,1)[4]", period = 4L,
    run = function() fit(log(data$UKgas), c(1L, 1L, 1L),
      list(order = c(1L, 1L, 1L), period = 4L))),
  conditional_seasonal_ar_nottem = list(
    dataset = "nottem", model = "AR(2) SAR(1)[12]", period = 12L,
    run = function() fit(data$nottem, c(2L, 0L, 0L),
      list(order = c(1L, 0L, 0L), period = 12L))),
  conditional_arma21_sunspots = list(
    dataset = "sunspots", model = "ARMA(2,1)", period = 12L,
    run = function() fit(data$sunspots, c(2L, 0L, 1L))),
  exact_ar2_nile = list(
    dataset = "Nile", model = "exact AR(2)", period = 1L,
    run = function() fit_ml(data$Nile, c(2L, 0L, 0L)))
)

if (length(selected)) cases <- cases[names(cases) %in% trimws(selected)]
if (!length(cases)) stop("no benchmark cases selected")

cat("runtime,benchmark,dataset,model,observations,period,warmups,measurements,repetitions,median_seconds_per_fit,converged\n")
for (name in names(cases)) {
  benchmark <- cases[[name]]
  for (iteration in seq_len(warmups)) benchmark$run()
  seconds <- numeric(measurements)
  result <- NULL
  for (iteration in seq_len(measurements)) {
    started <- proc.time()[["elapsed"]]
    for (repetition in seq_len(repetitions)) result <- benchmark$run()
    seconds[[iteration]] <- (proc.time()[["elapsed"]] - started) / repetitions
  }
  values <- data[[benchmark$dataset]]
  cat(paste("R", name, benchmark$dataset,
    encodeString(benchmark$model, quote = '"'), length(values), benchmark$period,
    warmups, measurements, repetitions, sprintf("%.6f", median(seconds)),
    result$code == 0L, sep = ","), "\n", sep = "")
}
