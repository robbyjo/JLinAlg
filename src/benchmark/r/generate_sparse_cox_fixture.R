#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

arguments <- commandArgs(trailingOnly = TRUE)
options <- list(
  output_dir = "build/benchmarks/sparse-cox-validation",
  individuals = 200L,
  observations_per_individual = 4L,
  genes = 8L,
  batches = 20L,
  seed = 20260905L
)
if (length(arguments) %% 2L) stop("arguments must be --name value pairs")
for (index in seq(1L, length(arguments), by = 2L))
  options[[sub("^--", "", arguments[[index]])]] <- arguments[[index + 1L]]
for (name in c("individuals", "observations_per_individual", "genes",
    "batches", "seed")) options[[name]] <- as.integer(options[[name]])
if (options$individuals < 40L || options$observations_per_individual < 1L ||
    options$genes < 1L || options$batches < 2L)
  stop("fixture dimensions are too small")

set.seed(options$seed)
dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
individual_ids <- sprintf("id%04d", seq_len(options$individuals))
founders <- min(40L, options$individuals)
sire <- rep(NA_character_, options$individuals)
dam <- rep(NA_character_, options$individuals)
if (options$individuals > founders) {
  descendants <- seq.int(founders + 1L, options$individuals)
  sire[descendants] <- individual_ids[1L + (descendants %% (founders / 2L))]
  dam[descendants] <- individual_ids[founders / 2L + 1L +
    ((descendants * 7L) %% (founders / 2L))]
}
pedigree <- data.table(id = individual_ids, sire = sire, dam = dam)

rows <- options$individuals * options$observations_per_individual
individual_index <- rep(seq_len(options$individuals),
  each = options$observations_per_individual)
observation_index <- sequence(rep(options$observations_per_individual,
  options$individuals))
batch_index <- 1L + ((individual_index * 11L + observation_index * 3L) %%
  options$batches)
analysis <- data.table(
  SampleName = sprintf("sample%06d", seq_len(rows)),
  Batch = sprintf("batch%03d", batch_index),
  animal_id = individual_ids[individual_index],
  Sex = as.integer(individual_index %% 2L),
  Age = 35 + ((individual_index * 13L + observation_index * 5L) %% 45L),
  BMI = 18 + ((individual_index * 17L + observation_index * 2L) %% 170L) / 10,
  WBC_Pred = sin(individual_index * 0.17) + observation_index * 0.03,
  LY_PER_Pred = cos(individual_index * 0.11) - observation_index * 0.02,
  MO_PER_Pred = sin(individual_index * 0.07 + observation_index),
  EO_PER_Pred = cos(individual_index * 0.05 + observation_index * 0.3),
  BA_PER_Pred = ((individual_index * 19L + observation_index) %% 101L) / 100
)
feature_keys <- sprintf("GENE%03d", seq_len(options$genes))
for (gene in seq_along(feature_keys)) {
  phase <- gene * 0.37
  analysis[, (feature_keys[[gene]]) :=
    sin(individual_index * (0.013 * gene + 0.021) + phase) +
    cos(seq_len(rows) * (0.007 * gene + 0.003)) +
    ((individual_index * (gene + 3L)) %% 17L) / 25]
}
features <- data.table(feature_key = feature_keys,
  feature_id = sprintf("synthetic-gene-%03d", seq_len(options$genes)))

standardized_age <- as.numeric(scale(analysis$Age))
standardized_bmi <- as.numeric(scale(analysis$BMI))
linear_predictor <- 0.25 * analysis$Sex + 0.32 * standardized_age -
  0.18 * standardized_bmi + 0.20 * sin(batch_index * 0.61)
event_time <- -log(runif(rows)) / exp(linear_predictor)
censor_time <- -log(runif(rows)) / 0.65
survival <- data.table(SampleName = analysis$SampleName,
  time = pmin(event_time, censor_time) + seq_len(rows) * 1e-12,
  event = as.integer(event_time <= censor_time))

fwrite(analysis, file.path(options$output_dir, "analysis.csv"))
fwrite(features, file.path(options$output_dir, "features.csv"))
fwrite(pedigree, file.path(options$output_dir, "pedigree.csv"), na = "")
fwrite(survival, file.path(options$output_dir, "survival.csv"))
manifest <- data.table(
  metric = c("seed", "individuals", "observations", "events", "genes",
    "batches", "observations_per_individual"),
  value = c(options$seed, options$individuals, rows, sum(survival$event),
    options$genes, options$batches, options$observations_per_individual)
)
fwrite(manifest, file.path(options$output_dir, "fixture_manifest.csv"))
print(manifest)
