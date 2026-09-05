# Deterministic comparator for the exact direct-surface JLinAlg LOESS path.
args <- commandArgs(trailingOnly = TRUE)
rows <- if (length(args) >= 1) as.integer(args[[1]]) else 5000L
span <- if (length(args) >= 2) as.numeric(args[[2]]) else 0.2
measurements <- if (length(args) >= 3) as.integer(args[[3]]) else 5L
warmups <- if (length(args) >= 4) as.integer(args[[4]]) else 2L
index <- seq_len(rows)
x <- -10 + 20 * (index - 1) / (rows - 1)
y <- sin(0.7 * x) + 0.05 * x + 0.1 * cos(0.37 * index)
control <- loess.control(surface = "direct", statistics = "none",
  trace.hat = "approximate")
fit_once <- function() loess(y ~ x, span = span, degree = 2,
  family = "gaussian", control = control, model = FALSE)
for (iteration in seq_len(warmups)) fit_once()
fit <- NULL
elapsed <- replicate(measurements,
  system.time(fit <<- fit_once())[["elapsed"]])
median_seconds <- median(elapsed)
checksum <- fit$fitted[[1]] + fit$fitted[[floor(rows / 2) + 1L]] +
  fit$fitted[[rows]]
cat("runtime,rows,span,degree,median_seconds,rows_per_second,checksum\n")
cat(sprintf("R_stats_direct,%d,%.6f,2,%.6f,%.2f,%.12f\n", rows, span,
  median_seconds, rows / median_seconds, checksum))
