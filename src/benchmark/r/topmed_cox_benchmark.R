#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(survival)
  library(coxme)
})

parse_arguments <- function(arguments) {
  result <- list(
    prepared_dir = "build/benchmarks/topmed-obesity",
    genes = 20L,
    max_rows = .Machine$integer.max,
    measurements = 3L,
    models = "cox,coxme,pedigree",
    fixed_variances = NULL,
    sparse_calc = NULL,
    pedigree_covariance = "sparse",
    output_prefix = NULL
  )
  index <- 1L
  while (index <= length(arguments)) {
    key <- sub("^--", "", arguments[[index]])
    if (index == length(arguments)) stop("arguments must be --name value pairs")
    result[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  result$genes <- as.integer(result$genes)
  result$max_rows <- as.integer(result$max_rows)
  result$measurements <- as.integer(result$measurements)
  if (!is.null(result$fixed_variances))
    result$fixed_variances <- as.numeric(strsplit(
      result$fixed_variances, ",", fixed = TRUE)[[1L]])
  if (!is.null(result$sparse_calc)) result$sparse_calc <-
    as.integer(result$sparse_calc)
  result
}

build_relationship <- function(pedigree, observed_ids) {
  ids <- as.character(pedigree$id)
  count <- length(ids)
  index <- setNames(seq_len(count), ids)
  sire <- match(as.character(pedigree$sire), ids)
  dam <- match(as.character(pedigree$dam), ids)
  sire[is.na(pedigree$sire) | pedigree$sire == ""] <- NA_integer_
  dam[is.na(pedigree$dam) | pedigree$dam == ""] <- NA_integer_

  parent <- seq_len(count)
  find_root <- function(value) {
    while (parent[[value]] != value) {
      parent[[value]] <<- parent[[parent[[value]]]]
      value <- parent[[value]]
    }
    value
  }
  unite <- function(left, right) {
    first <- find_root(left)
    second <- find_root(right)
    if (first != second) parent[[second]] <<- first
  }
  children <- vector("list", count)
  indegree <- integer(count)
  for (individual in seq_len(count)) {
    for (relative in c(sire[[individual]], dam[[individual]])) {
      if (!is.na(relative)) {
        unite(individual, relative)
        children[[relative]] <- c(children[[relative]], individual)
        indegree[[individual]] <- indegree[[individual]] + 1L
      }
    }
  }
  queue <- which(indegree == 0L)
  order <- integer(count)
  head <- 1L
  tail <- length(queue)
  emitted <- 0L
  while (head <= tail) {
    value <- queue[[head]]
    head <- head + 1L
    emitted <- emitted + 1L
    order[[emitted]] <- value
    for (child in children[[value]]) {
      indegree[[child]] <- indegree[[child]] - 1L
      if (indegree[[child]] == 0L) {
        tail <- tail + 1L
        queue[[tail]] <- child
      }
    }
  }
  if (emitted != count) stop("pedigree contains a cycle")
  roots <- vapply(seq_len(count), find_root, integer(1L))
  components <- split(order, roots[order])

  matrix_i <- integer(0L)
  matrix_j <- integer(0L)
  matrix_x <- numeric(0L)
  for (component in components) {
    size <- length(component)
    local <- setNames(seq_len(size), component)
    relationship <- matrix(0, size, size)
    for (position in seq_len(size)) {
      individual <- component[[position]]
      sire_position <- if (is.na(sire[[individual]])) NA_integer_
        else unname(local[[as.character(sire[[individual]])]])
      dam_position <- if (is.na(dam[[individual]])) NA_integer_
        else unname(local[[as.character(dam[[individual]])]])
      if (!is.na(sire_position) && !is.na(dam_position)) {
        if (position > 1L) relationship[position, seq_len(position - 1L)] <-
          0.5 * (relationship[sire_position, seq_len(position - 1L)] +
            relationship[dam_position, seq_len(position - 1L)])
        relationship[position, position] <-
          1 + 0.5 * relationship[sire_position, dam_position]
      } else if (!is.na(sire_position) || !is.na(dam_position)) {
        known <- if (!is.na(sire_position)) sire_position else dam_position
        if (position > 1L) relationship[position, seq_len(position - 1L)] <-
          0.5 * relationship[known, seq_len(position - 1L)]
        relationship[position, position] <- 1
      } else {
        relationship[position, position] <- 1
      }
      if (position > 1L)
        relationship[seq_len(position - 1L), position] <-
          relationship[position, seq_len(position - 1L)]
    }
    entries <- which(relationship != 0, arr.ind = TRUE)
    matrix_i <- c(matrix_i, component[entries[, 1L]])
    matrix_j <- c(matrix_j, component[entries[, 2L]])
    matrix_x <- c(matrix_x, relationship[entries])
  }
  result <- sparseMatrix(i = matrix_i, j = matrix_j, x = matrix_x,
    dims = c(count, count), dimnames = list(ids, ids), symmetric = FALSE)
  selected <- unique(as.character(observed_ids))
  forceSymmetric(result[selected, selected, drop = FALSE], uplo = "L")
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
setDTthreads(1L)
cat(sprintf("R=%s survival=%s coxme=%s data.table=%s\n",
  paste(R.version$major, R.version$minor, sep = "."),
  as.character(packageVersion("survival")),
  as.character(packageVersion("coxme")),
  as.character(packageVersion("data.table"))))
analysis <- fread(file.path(options$prepared_dir, "analysis.csv"))
survival_path <- file.path(options$prepared_dir, "survival.csv")
if (!file.exists(survival_path)) {
  set.seed(20260903)
  standardized_age <- as.numeric(scale(analysis$Age))
  standardized_bmi <- as.numeric(scale(analysis$BMI))
  levy_index <- as.integer(factor(analysis$Batch))
  linear_predictor <- 0.25 * analysis$Sex + 0.35 * standardized_age +
    0.20 * standardized_bmi + 0.30 * sin(1.7 * levy_index)
  event_time <- -log(runif(nrow(analysis))) / exp(linear_predictor)
  censor_time <- -log(runif(nrow(analysis))) / 0.7
  generated <- data.table(SampleName = analysis$SampleName,
    time = pmin(event_time, censor_time) + 1e-8,
    event = as.integer(event_time <= censor_time))
  fwrite(generated, survival_path)
}
if (nrow(analysis) > options$max_rows) {
  selected_rows <- unique(as.integer(round(seq(
    1, nrow(analysis), length.out = options$max_rows))))
  analysis <- analysis[selected_rows]
}
survival_data <- fread(survival_path)
features <- fread(file.path(options$prepared_dir, "features.csv"),
  nrows = options$genes)
if (nrow(features) != options$genes) stop("not enough prepared genes")
analysis <- merge(analysis, survival_data, by = "SampleName", sort = FALSE)
analysis[, Levy_Set := factor(Batch)]
analysis[, animal_id := factor(animal_id)]
variable <- vapply(features$feature_key,
  function(gene) length(unique(analysis[[gene]])) > 1L, logical(1L))
if (any(!variable))
  cat(sprintf("excluded constant genes=%d\n", sum(!variable)))
features <- features[variable]
cat(sprintf("R cohort rows=%d events=%d genes=%d\n",
  nrow(analysis), sum(analysis$event), nrow(features)))
requested_models <- strsplit(options$models, ",", fixed = TRUE)[[1L]]
coxme_control <- if (is.null(options$sparse_calc)) coxme.control() else
  coxme.control(sparse.calc = options$sparse_calc)
cat(sprintf("R coxme sparse.calc=%s\n",
  if (is.null(options$sparse_calc)) "automatic" else options$sparse_calc))

fixed_terms <- c("GENE", "Sex", "Age", "WBC_Pred", "LY_PER_Pred",
  "MO_PER_Pred", "EO_PER_Pred", "BA_PER_Pred")
fixed_formula <- function(gene) as.formula(paste(
  "Surv(time, event) ~", paste(sub("GENE", gene, fixed_terms), collapse = " + ")))
mixed_formula <- function(gene) update(fixed_formula(gene), . ~ . + (1 | Levy_Set))
pedigree_formula <- function(gene) update(fixed_formula(gene),
  . ~ . + (1 | animal_id) + (1 | Levy_Set))

pedigree_variance <- NULL
if (any(c("pedigree", "pedigree-dense") %in% requested_models)) {
  pedigree <- fread(file.path(options$prepared_dir, "pedigree.csv"),
    colClasses = "character", na.strings = c("", "NA"))
  pedigree_variance <- build_relationship(pedigree, levels(analysis$animal_id))
  cat(sprintf("R pedigree covariance dimension=%d nonzeros=%d\n",
    nrow(pedigree_variance), length(pedigree_variance@x)))
  if (identical(Sys.getenv("JLINALG_COX_VALIDATE_PRECISION"), "true") &&
      setequal(as.character(pedigree$id), levels(analysis$animal_id))) {
    ids <- as.character(pedigree$id)
    index <- setNames(seq_along(ids), ids)
    precision <- matrix(0, length(ids), length(ids))
    for (individual in seq_along(ids)) {
      sire <- if (is.na(pedigree$sire[[individual]])) NA_integer_ else
        unname(index[[as.character(pedigree$sire[[individual]])]])
      dam <- if (is.na(pedigree$dam[[individual]])) NA_integer_ else
        unname(index[[as.character(pedigree$dam[[individual]])]])
      mendelian <- if (is.na(sire) && is.na(dam)) 1 else
        if (is.na(sire) || is.na(dam)) 0.75 else 0.5
      relatives <- c(individual, sire, dam)
      coefficients <- c(1, -0.5, -0.5)
      retained <- !is.na(relatives)
      precision[relatives[retained], relatives[retained]] <-
        precision[relatives[retained], relatives[retained]] +
        tcrossprod(coefficients[retained]) / mendelian
    }
    covariance_ordered <- as.matrix(pedigree_variance[ids, ids])
    cat(sprintf("R pedigree precision parity max_abs=%.12g\n",
      max(abs(solve(covariance_ordered) - precision))))
  }
  if (options$pedigree_covariance == "dense") {
    pedigree_variance <- as.matrix(pedigree_variance)
  } else if (options$pedigree_covariance != "sparse") {
    stop("--pedigree_covariance must be sparse or dense")
  }
  cat(sprintf("R pedigree covariance representation=%s\n",
    options$pedigree_covariance))
}

fit_gene <- function(model, gene) {
  if (model == "cox") {
    coxph(fixed_formula(gene), data = analysis, ties = "efron",
      model = FALSE, x = FALSE, y = FALSE)
  } else if (model == "coxme") {
    coxme(mixed_formula(gene), data = analysis, ties = "efron",
      control = coxme_control,
      vfixed = if (is.null(options$fixed_variances)) NULL
        else options$fixed_variances[[length(options$fixed_variances)]])
  } else if (model %in% c("pedigree", "pedigree-dense")) {
    coxme(pedigree_formula(gene), data = analysis, ties = "efron",
      control = coxme_control,
      varlist = list(coxmeMlist(pedigree_variance, rescale = FALSE),
        coxmeFull()), vfixed = if (is.null(options$fixed_variances)) NULL
          else as.list(options$fixed_variances))
  } else stop(paste("unknown model", model))
}

timings <- list()
results <- list()
genes <- features$feature_key
for (model in requested_models) {
  warm_fits <- lapply(genes, function(gene) fit_gene(model, gene))
  if (model != "cox")
    cat(sprintf("R model=%s variance_terms=%s\n", model,
      paste(names(warm_fits[[1L]]$vcoef), collapse = ";")))
  for (measurement in seq_len(options$measurements)) {
    gc()
    elapsed <- system.time(fits <- lapply(genes, function(gene)
      fit_gene(model, gene)))[["elapsed"]]
    timing <- data.table(runtime = "R", model = model, threads = 1L,
      measurement = measurement, genes = length(genes), rows = nrow(analysis),
      seconds = elapsed,
      genes_per_second = length(genes) / elapsed)
    timings[[length(timings) + 1L]] <- timing
    print(timing)
    if (measurement == 1L) {
      for (gene_index in seq_along(genes)) {
        fit <- fits[[gene_index]]
        beta <- if (model == "cox") coef(fit) else fixef(fit)
        standard_error <- sqrt(diag(vcov(fit)))
        coefficient <- match(genes[[gene_index]], names(beta))
        frailty_variances <- if (model == "cox") "" else
          paste(as.numeric(unlist(fit$vcoef)), collapse = ";")
        results[[length(results) + 1L]] <- data.table(
          runtime = "R", model = model, threads = 1L,
          feature_key = genes[[gene_index]],
          feature_id = features$feature_id[[gene_index]],
          beta = beta[[coefficient]],
          standard_error = standard_error[[coefficient]],
          frailty_variances = frailty_variances)
      }
    }
  }
}

output_prefix <- if (is.null(options$output_prefix))
  file.path(options$prepared_dir, "r_cox") else options$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(timings), paste0(output_prefix, "_timings.csv"))
fwrite(rbindlist(results), paste0(output_prefix, "_results.csv"))
