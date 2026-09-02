# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

suppressPackageStartupMessages(library(data.table))

parse_arguments <- function(arguments) {
  result <- list(
    mode = "benchmark",
    data_dir = "D:/Research/topmed/splicing-bmi/new",
    prepared_dir = "build/benchmarks/topmed100",
    genes = 100L,
    measurements = 1L,
    models = "ols,reml,pedigree"
  )
  index <- 1L
  while (index <= length(arguments)) {
    key <- sub("^--", "", arguments[[index]])
    if (index == length(arguments)) stop("missing value for --", key)
    value <- arguments[[index + 1L]]
    key <- gsub("-", "_", key, fixed = TRUE)
    result[[key]] <- value
    index <- index + 2L
  }
  result$genes <- as.integer(result$genes)
  result$measurements <- as.integer(result$measurements)
  result
}

require_column <- function(table, requested) {
  index <- match(tolower(requested), tolower(names(table)))
  if (is.na(index)) stop("missing column: ", requested)
  names(table)[[index]]
}

as_identifier <- function(value) {
  result <- trimws(as.character(value))
  result[result == "" | toupper(result) %in% c("NA", "N/A", ".")] <- NA_character_
  result
}

prepare_inputs <- function(options) {
  dir.create(options$prepared_dir, recursive = TRUE, showWarnings = FALSE)
  phenotype_path <- file.path(options$data_dir, "mastermat-batch1234-wcbc-ext.csv")
  omics_path <- file.path(options$data_dir, "fhs-batch1234-jointtmmfpkm-genes-log2.csv")
  pedigree_path <- file.path(options$data_dir, "sabre_ped_0407_v1_rj.csv")
  annotation_path <- file.path(options$data_dir, "gencode.v48.annotation.gene.csv")

  phenotype <- fread(phenotype_path,
    colClasses = list(character = c("SampleName", "sabreid", "Levy_Set")),
    showProgress = FALSE)
  sample_column <- require_column(phenotype, "SAMPLENAME")
  requested <- c("BMI", "Obesity", "Sex", "Age", "WBC_Pred", "LY_PER_Pred",
    "MO_PER_Pred", "EO_PER_Pred", "BA_PER_Pred", "Levy_Set", "sabreid")
  actual <- vapply(requested, function(value) require_column(phenotype, value), "")
  setnames(phenotype, sample_column, "SampleName")
  actual[actual == sample_column] <- "SampleName"
  phenotype[, SampleName := as_identifier(SampleName)]
  phenotype[, sabreid := as_identifier(get(actual[["sabreid"]]))]

  header <- strsplit(readLines(omics_path, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1L]]
  omics_ids <- header[-1L]
  phenotype_index <- match(omics_ids, phenotype$SampleName)
  matched <- !is.na(phenotype_index)
  complete <- rep(FALSE, length(omics_ids))
  complete[matched] <- complete.cases(phenotype[phenotype_index[matched], ..actual])
  retained_omics_columns <- which(complete) + 1L
  retained_phenotype <- phenotype[phenotype_index[complete]]
  if (nrow(retained_phenotype) < 10L) stop("too few complete matched samples")

  candidate_count <- max(500L, options$genes * 3L)
  expression <- fread(omics_path, skip = 1L, header = FALSE,
    nrows = candidate_count, select = c(1L, retained_omics_columns),
    showProgress = FALSE)
  setnames(expression, c("feature_id", retained_phenotype$SampleName))
  feature_id <- expression[["feature_id"]]
  values <- as.matrix(expression[, -1L])
  storage.mode(values) <- "double"
  valid <- apply(values, 1L, function(row) all(is.finite(row)) && sd(row) > 0)
  selected <- which(valid)[seq_len(min(sum(valid), options$genes))]
  if (length(selected) != options$genes) {
    stop("only ", length(selected), " complete nonconstant genes among first ",
      candidate_count, " candidates")
  }
  values <- values[selected, , drop = FALSE]
  feature_id <- feature_id[selected]
  feature_key <- sprintf("gene_%03d", seq_along(feature_id))

  annotation <- fread(annotation_path, select = c("gene_id", "gene_name"),
    colClasses = "character", showProgress = FALSE)
  annotation <- unique(annotation, by = "gene_id")
  gene_name <- annotation$gene_name[match(feature_id, annotation$gene_id)]
  metadata <- data.table(feature_key, feature_id, gene_name)

  analysis <- data.table(
    SampleName = retained_phenotype$SampleName,
    animal_id = retained_phenotype$sabreid,
    Batch = as.character(retained_phenotype[[actual[["Levy_Set"]]]]),
    BMI = as.numeric(retained_phenotype[[actual[["BMI"]]]]),
    Obesity = as.numeric(retained_phenotype[[actual[["Obesity"]]]]),
    Sex = as.numeric(retained_phenotype[[actual[["Sex"]]]]),
    Age = as.numeric(retained_phenotype[[actual[["Age"]]]]),
    WBC_Pred = as.numeric(retained_phenotype[[actual[["WBC_Pred"]]]]),
    LY_PER_Pred = as.numeric(retained_phenotype[[actual[["LY_PER_Pred"]]]]),
    MO_PER_Pred = as.numeric(retained_phenotype[[actual[["MO_PER_Pred"]]]]),
    EO_PER_Pred = as.numeric(retained_phenotype[[actual[["EO_PER_Pred"]]]]),
    BA_PER_Pred = as.numeric(retained_phenotype[[actual[["BA_PER_Pred"]]]])
  )
  for (column in seq_along(feature_key)) {
    set(analysis, j = feature_key[[column]], value = values[column, ])
  }

  pedigree <- fread(pedigree_path,
    colClasses = list(character = c("sabreid", "fid", "mid")),
    showProgress = FALSE)
  id_column <- require_column(pedigree, "sabreid")
  sire_column <- require_column(pedigree, "fid")
  dam_column <- require_column(pedigree, "mid")
  pedigree <- pedigree[, .(
    id = as_identifier(get(id_column)),
    sire = as_identifier(get(sire_column)),
    dam = as_identifier(get(dam_column)))]
  pedigree <- pedigree[!is.na(id)]
  if (anyDuplicated(pedigree$id)) stop("pedigree IDs are not unique")

  original_ids <- pedigree$id
  selected_ids <- intersect(unique(analysis$animal_id), original_ids)
  frontier <- selected_ids
  while (length(frontier)) {
    rows <- match(frontier, pedigree$id)
    parents <- unique(c(pedigree$sire[rows], pedigree$dam[rows]))
    parents <- parents[!is.na(parents) & !(parents %in% selected_ids)]
    selected_ids <- c(selected_ids, parents)
    frontier <- parents
  }
  pedigree <- pedigree[id %in% selected_ids]
  singleton_ids <- setdiff(unique(analysis$animal_id), original_ids)
  if (length(singleton_ids)) {
    pedigree <- rbind(pedigree,
      data.table(id = singleton_ids, sire = NA_character_, dam = NA_character_))
  }

  ordered_ids <- character()
  remaining <- pedigree$id
  while (length(remaining)) {
    rows <- match(remaining, pedigree$id)
    ready <- (is.na(pedigree$sire[rows]) | pedigree$sire[rows] %in% ordered_ids) &
      (is.na(pedigree$dam[rows]) | pedigree$dam[rows] %in% ordered_ids)
    if (!any(ready)) stop("pedigree ancestry cycle or missing parent")
    ordered_ids <- c(ordered_ids, remaining[ready])
    remaining <- remaining[!ready]
  }
  pedigree <- pedigree[match(ordered_ids, id)]

  fwrite(analysis, file.path(options$prepared_dir, "analysis.csv"))
  fwrite(metadata, file.path(options$prepared_dir, "features.csv"))
  fwrite(pedigree, file.path(options$prepared_dir, "pedigree.csv"), na = "")
  manifest <- data.table(
    metric = c("samples", "genes", "pedigree_individuals", "singletons",
      "batch_source", "sample_id_source", "pedigree_id_source"),
    value = c(nrow(analysis), nrow(metadata), nrow(pedigree), length(singleton_ids),
      "Levy_Set", "SampleName (case-insensitive SAMPLENAME match)", "sabreid"))
  fwrite(manifest, file.path(options$prepared_dir, "manifest.csv"))
  print(manifest)
}

lmer_control <- function() {
  lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE,
    check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "ignore",
    check.nobs.vs.nRE = "ignore", check.rankX = "silent.drop.cols",
    check.conv.singular = lme4:::.makeCC(action = "ignore", tol = 1e-4))
}

benchmark_models <- function(options) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(pedigreemm))
  analysis <- fread(file.path(options$prepared_dir, "analysis.csv"),
    colClasses = list(character = c("SampleName", "animal_id", "Batch")),
    showProgress = FALSE)
  features <- fread(file.path(options$prepared_dir, "features.csv"),
    colClasses = "character", showProgress = FALSE)
  pedigree_table <- fread(file.path(options$prepared_dir, "pedigree.csv"),
    colClasses = "character", na.strings = "", showProgress = FALSE)
  features <- head(features, options$genes)
  analysis[, Batch := factor(Batch)]
  analysis[, animal_id := factor(animal_id, levels = pedigree_table$id)]
  pedigree <- pedigreemm::pedigree(sire = pedigree_table$sire,
    dam = pedigree_table$dam, label = pedigree_table$id)
  pedigree_list <- list(animal_id = pedigree)
  controls <- lmer_control()
  fixed_rhs <- "Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred"
  requested_models <- strsplit(options$models, ",", fixed = TRUE)[[1L]]

  fit_gene <- function(model, feature) {
    formula <- as.formula(paste("BMI ~", feature, "+", fixed_rhs,
      if (model == "reml") "+ (1|Batch)" else
      if (model == "pedigree") "+ (1|Batch) + (1|animal_id)" else ""))
    if (model == "ols") return(lm(formula, data = analysis))
    if (model == "reml") return(lmer(formula, data = analysis,
      REML = TRUE, control = controls))
    if (model == "pedigree") return(pedigreemm(formula, data = analysis,
      pedigree = pedigree_list, control = controls))
    stop("unknown model: ", model)
  }

  timing_rows <- list()
  result_rows <- list()
  for (model in requested_models) {
    invisible(fit_gene(model, features$feature_key[[1L]]))
    for (measurement in seq_len(options$measurements)) {
      garbage <- gc()
      elapsed <- system.time({
        fits <- lapply(features$feature_key, function(feature) fit_gene(model, feature))
      })[["elapsed"]]
      timing_rows[[length(timing_rows) + 1L]] <- data.table(
        runtime = "R", model, threads = 1L, measurement,
        genes = nrow(features), seconds = elapsed,
        genes_per_second = nrow(features) / elapsed)
      if (measurement == 1L) {
        for (index in seq_along(fits)) {
          coefficients <- if (model == "ols") coef(summary(fits[[index]])) else
            coef(summary(fits[[index]]))
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
  fwrite(timings, file.path(options$prepared_dir, "r_timings.csv"))
  fwrite(results, file.path(options$prepared_dir, "r_results.csv"))
  print(timings)
  if (length(warnings())) print(warnings())
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (options$mode == "prepare") {
  prepare_inputs(options)
} else if (options$mode == "benchmark") {
  benchmark_models(options)
} else {
  stop("--mode must be prepare or benchmark")
}
