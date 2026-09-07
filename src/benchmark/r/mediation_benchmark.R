# Deterministic Gaussian mediation comparator for JLinAlg.
# The non-sampling path is three lm() fits; lme4 is optional for the sparse
# mixed-effect comparison. Run with: Rscript mediation_benchmark.R [rows]
args <- commandArgs(trailingOnly = TRUE)
rows <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
groups <- if (length(args) >= 2) as.integer(args[[2]]) else 100L
measurements <- if (length(args) >= 3) as.integer(args[[3]]) else 3L
warmups <- if (length(args) >= 4) as.integer(args[[4]]) else 1L
index <- seq_len(rows)
treatment <- -1 + 2 * (index - 1) / (rows - 1)
covariate <- sin(0.17 * index)
group <- factor((index - 1) %% groups)
group_effect <- 0.4 * sin(0.7 * ((index - 1) %% groups))
mediator <- 0.8 + 0.6 * treatment + 0.4 * covariate + group_effect +
  0.1 * cos(0.31 * index)
outcome <- 0.3 + 0.25 * treatment + 1.1 * mediator + 0.2 * covariate +
  group_effect + 0.1 * sin(0.23 * index)
data <- data.frame(outcome, treatment, mediator, covariate, group)

fit_ols <- function() {
  mediator_model <- lm(mediator ~ treatment + covariate, data = data)
  outcome_model <- lm(outcome ~ treatment + mediator + covariate, data = data)
  total_model <- lm(outcome ~ treatment + covariate, data = data)
  c(a = coef(mediator_model)[["treatment"]],
    b = coef(outcome_model)[["mediator"]],
    indirect = coef(mediator_model)[["treatment"]] *
      coef(outcome_model)[["mediator"]],
    direct = coef(outcome_model)[["treatment"]],
    total = coef(total_model)[["treatment"]])
}

measure <- function(name, fit_once) {
  for (iteration in seq_len(warmups)) fit_once()
  elapsed <- numeric(measurements)
  result <- NULL
  for (iteration in seq_len(measurements)) {
    started <- proc.time()[["elapsed"]]
    result <- fit_once()
    elapsed[[iteration]] <- proc.time()[["elapsed"]] - started
  }
  median_seconds <- median(elapsed)
  cat(sprintf("R,%s,%d,%d,%.6f,%.2f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",
    name, rows, groups, median_seconds, rows / median_seconds,
    result[["a"]], result[["b"]], result[["indirect"]],
    result[["direct"]], result[["total"]], sum(result)))
}

cat("runtime,model,rows,groups,median_seconds,rows_per_second,a,b,indirect,direct,total,checksum\n")
measure("ols", fit_ols)

if (requireNamespace("lme4", quietly = TRUE)) {
  fit_mixed <- function() {
    mediator_model <- lme4::lmer(
      mediator ~ treatment + covariate + (1 | group), data = data, REML = TRUE)
    outcome_model <- lme4::lmer(
      outcome ~ treatment + mediator + covariate + (1 | group),
      data = data, REML = TRUE)
    total_model <- lme4::lmer(
      outcome ~ treatment + covariate + (1 | group), data = data, REML = TRUE)
    c(a = lme4::fixef(mediator_model)[["treatment"]],
      b = lme4::fixef(outcome_model)[["mediator"]],
      indirect = lme4::fixef(mediator_model)[["treatment"]] *
        lme4::fixef(outcome_model)[["mediator"]],
      direct = lme4::fixef(outcome_model)[["treatment"]],
      total = lme4::fixef(total_model)[["treatment"]])
  }
  measure("mixed_lmer", fit_mixed)
} else {
  message("lme4 is unavailable; skipped mixed_lmer comparison")
}
