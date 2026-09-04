args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[[1L]] else
  "D:/Research/topmed/splicing-bmi/new/mastermat-batch1234-wcbc-ext.csv"
measurements <- if (length(args) >= 2L) as.integer(args[[2L]]) else 10L
local_library <- normalizePath(".r-reference-lib", winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))
suppressPackageStartupMessages(library(lavaan))

variables <- c(
  "Sex", "Age", "BMI", "Waist", "Systolic_BP", "Diastolic_BP",
  "Glucose", "HDL", "LnTG", "LogInsulin", "CRP", "eGFR"
)
data <- read.csv(input, check.names = FALSE)[variables]
data <- data[complete.cases(data), , drop = FALSE]
data[] <- lapply(data, function(value) as.numeric(scale(value)))

model <- '
  BMI ~ Age + Sex
  Waist ~ BMI + Sex
  Systolic_BP ~ Age + BMI
  Diastolic_BP ~ Age + BMI
  Glucose ~ Age + BMI
  HDL ~ BMI + LnTG + Sex
  LnTG ~ BMI + Glucose + Sex
  LogInsulin ~ BMI + Glucose
  CRP ~ BMI + Sex
  eGFR ~ Age + Sex + BMI
  Age ~~ Sex
'
fit_once <- function() sem(
  model, data = data, meanstructure = FALSE, fixed.x = FALSE,
  auto.cov.y = FALSE, estimator = "ML", information = "expected",
  se = "standard"
)
invisible(fit_once())
times <- numeric(measurements)
fit <- NULL
for (measurement in seq_len(measurements)) {
  started <- proc.time()[["elapsed"]]
  fit <- fit_once()
  times[[measurement]] <- proc.time()[["elapsed"]] - started
}
measures <- fitMeasures(
  fit, c("logl", "chisq", "df", "cfi", "tli", "rmsea", "srmr"))
parameters <- parameterEstimates(fit)
parameters <- parameters[parameters$op %in% c("~", "~~") &
  !(parameters$op == "~~" & parameters$lhs != parameters$rhs &
    !((parameters$lhs == "Age" & parameters$rhs == "Sex") |
      (parameters$lhs == "Sex" & parameters$rhs == "Age"))), ]
labels <- ifelse(parameters$op == "~",
  paste0(parameters$lhs, "~", parameters$rhs),
  paste0(parameters$lhs, "~~", parameters$rhs))

cat("runtime,observations,variables,parameters,iterations,converged,",
    "median_seconds,log_likelihood,chi_square,df,cfi,tli,rmsea,srmr\n", sep = "")
cat(sprintf(
  "lavaan-%s,%d,%d,%d,%d,%s,%.9f,%.12g,%.12g,%d,%.12g,%.12g,%.12g,%.12g\n",
  as.character(packageVersion("lavaan")), nobs(fit), length(variables),
  nrow(parameters), lavInspect(fit, "iterations"),
  lavInspect(fit, "converged"), median(times), measures[["logl"]],
  measures[["chisq"]], as.integer(measures[["df"]]), measures[["cfi"]],
  measures[["tli"]], measures[["rmsea"]], measures[["srmr"]]
))
cat("label,estimate,se\n")
for (index in seq_len(nrow(parameters)))
  cat(sprintf("%s,%.12g,%.12g\n", labels[[index]],
    parameters$est[[index]], parameters$se[[index]]))
