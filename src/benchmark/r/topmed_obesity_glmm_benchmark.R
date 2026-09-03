# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

suppressPackageStartupMessages(library(data.table))

parse_arguments <- function(arguments) {
  result <- list(
    prepared_dir = "build/benchmarks/topmed100",
    output_prefix = NULL,
    genes = 20L,
    measurements = 3L,
    models = "glm,glmm,pedigree"
  )
  index <- 1L
  while (index <= length(arguments)) {
    key <- gsub("-", "_", sub("^--", "", arguments[[index]]), fixed = TRUE)
    if (index == length(arguments)) stop("missing value for --", key)
    result[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  result$genes <- as.integer(result$genes)
  result$measurements <- as.integer(result$measurements)
  result
}

glmer_control <- function() {
  lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE,
    check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "ignore",
    check.nobs.vs.nRE = "ignore", check.rankX = "silent.drop.cols",
    check.conv.singular = lme4:::.makeCC(action = "ignore", tol = 1e-4))
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(pedigreemm))

analysis <- fread(file.path(options$prepared_dir, "analysis.csv"),
  colClasses = list(character = c("SampleName", "animal_id", "Batch")),
  showProgress = FALSE)
features <- head(fread(file.path(options$prepared_dir, "features.csv"),
  colClasses = "character", showProgress = FALSE), options$genes)
pedigree_table <- fread(file.path(options$prepared_dir, "pedigree.csv"),
  colClasses = "character", na.strings = "", showProgress = FALSE)
analysis[, Batch := factor(Batch)]
analysis[, animal_id := factor(animal_id, levels = pedigree_table$id)]
pedigree <- pedigreemm::pedigree(sire = pedigree_table$sire,
  dam = pedigree_table$dam, label = pedigree_table$id)
pedigree_list <- list(animal_id = pedigree)
controls <- glmer_control()
fixed_rhs <- paste("Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred +",
  "EO_PER_Pred + BA_PER_Pred")
requested_models <- strsplit(options$models, ",", fixed = TRUE)[[1L]]

fit_gene <- function(model, feature) {
  random_rhs <- if (model == "glmm") "+ (1|Batch)" else
    if (model == "pedigree") "+ (1|Batch) + (1|animal_id)" else ""
  formula <- as.formula(paste("Obesity ~", feature, "+", fixed_rhs, random_rhs))
  if (model == "glm") {
    return(glm(formula, data = analysis, family = binomial()))
  }
  if (model == "glmm") {
    return(glmer(formula, data = analysis, family = binomial(), nAGQ = 1L,
      control = controls))
  }
  if (model == "pedigree") {
    return(pedigreemm(formula, data = analysis, family = binomial(),
      pedigree = pedigree_list, control = controls))
  }
  stop("unknown model: ", model)
}

timing_rows <- list()
result_rows <- list()
for (model in requested_models) {
  invisible(fit_gene(model, features$feature_key[[1L]]))
  for (measurement in seq_len(options$measurements)) {
    invisible(gc())
    elapsed <- system.time({
      fits <- lapply(features$feature_key,
        function(feature) fit_gene(model, feature))
    })[["elapsed"]]
    timing_rows[[length(timing_rows) + 1L]] <- data.table(
      runtime = "R", model, threads = 1L, measurement,
      genes = nrow(features), seconds = elapsed,
      genes_per_second = nrow(features) / elapsed)
    print(timing_rows[[length(timing_rows)]])
    if (measurement == 1L) {
      for (index in seq_along(fits)) {
        coefficients <- coef(summary(fits[[index]]))
        key <- features$feature_key[[index]]
        result_rows[[length(result_rows) + 1L]] <- data.table(
          runtime = "R", model, threads = 1L, feature_key = key,
          feature_id = features$feature_id[[index]],
          beta = coefficients[key, "Estimate"],
          standard_error = coefficients[key, "Std. Error"])
      }
    }
  }
}

timings <- rbindlist(timing_rows)
results <- rbindlist(result_rows)
output_prefix <- if (is.null(options$output_prefix))
  file.path(options$prepared_dir, "r_obesity_glmm") else options$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(timings, paste0(output_prefix, "_timings.csv"))
fwrite(results, paste0(output_prefix, "_results.csv"))
print(timings)
if (length(warnings())) print(warnings())
