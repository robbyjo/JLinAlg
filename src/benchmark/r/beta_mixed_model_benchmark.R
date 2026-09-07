# Comparable grouped random-intercept workload for glmmTMB.
.libPaths(c("build/r-library", .libPaths()))
suppressPackageStartupMessages(library(glmmTMB))

rows <- as.integer(Sys.getenv("BETA_MIXED_ROWS", "30000"))
groups_n <- as.integer(Sys.getenv("BETA_MIXED_GROUPS", "1000"))
repetitions <- as.integer(Sys.getenv("BETA_MIXED_REPETITIONS", "3"))
index <- 0:(rows - 1L)
group_index <- index %% groups_n
x <- -1 + 2 * ((index * 37) %% 1009) / 1008
random <- 0.55 * sin(1.7 * group_index) + 0.15 * cos(0.31 * group_index)
mu <- plogis(-0.3 + 0.8 * x + random)
standard_deviation <- sqrt(mu * (1 - mu) / 19)
residual <- (sin(1.73 * index + 0.2) + 0.55 * cos(0.47 * index)) / 1.14
y <- pmax(1e-5, pmin(1 - 1e-5, mu + standard_deviation * residual))
data <- data.frame(y = y, x = x, group = factor(group_index))

fit_once <- function() glmmTMB(y ~ x + (1 | group), data = data,
  family = beta_family(link = "logit"))
fit_once()
times <- numeric(repetitions)
fit <- NULL
for (iteration in seq_len(repetitions)) {
  start <- proc.time()[["elapsed"]]
  fit <- fit_once()
  times[iteration] <- (proc.time()[["elapsed"]] - start) * 1000
}
cat(sprintf("rows=%d groups=%d repetitions=%d\n", rows, groups_n, repetitions))
cat(sprintf(paste0("glmmTMB median_ms=%.3f beta=%s precision=%.6g ",
  "variance=%.6g\n"), median(times),
  paste(format(unname(fixef(fit)$cond), digits = 8), collapse = ","),
  exp(unname(fixef(fit)$disp[1])), unname(VarCorr(fit)$cond$group[1])))
