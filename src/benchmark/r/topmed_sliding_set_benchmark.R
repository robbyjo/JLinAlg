suppressPackageStartupMessages({
  local_lib <- normalizePath(".r-reference-lib", mustWork = FALSE)
  .libPaths(c(local_lib, .libPaths()))
  library(data.table)
  library(Matrix)
  library(GMMAT)
})

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")

options <- list(
  prepared_dir = "build/benchmarks/topmed-sliding-set",
  measurements = 3L,
  models = c("BMI", "Obesity"))
args <- commandArgs(trailingOnly = TRUE)
index <- 1L
while (index <= length(args)) {
  key <- args[[index]]
  if (index == length(args)) stop("missing value for ", key)
  value <- args[[index + 1L]]
  if (key == "--prepared-dir") options$prepared_dir <- value
  else if (key == "--measurements") options$measurements <- as.integer(value)
  else if (key == "--models") options$models <- strsplit(value, ",", fixed = TRUE)[[1L]]
  else stop("unknown option: ", key)
  index <- index + 2L
}

analysis <- fread(file.path(options$prepared_dir, "analysis.csv"),
  colClasses = list(character = c("framid", "sabreid")))
variants <- fread(file.path(options$prepared_dir, "variants.csv"))
windows <- fread(file.path(options$prepared_dir, "windows.csv"))
genotypes <- fread(file.path(options$prepared_dir, "genotypes.csv"),
  colClasses = list(character = "framid"))
stopifnot(identical(analysis$framid, genotypes$framid))
stopifnot(identical(names(genotypes)[-1L], variants$key))

n <- nrow(analysis)
connection <- file(file.path(options$prepared_dir, "relationship.bin"), "rb")
relationship <- matrix(readBin(connection, what = "double", n = n * n,
  size = 8L, endian = "big"), nrow = n, byrow = TRUE)
close(connection)
if (length(relationship) != n * n) stop("incomplete relationship matrix")
rownames(relationship) <- colnames(relationship) <- analysis$framid

# Orientation, imputation, weights, window membership, and null fitting are all
# prepared before the stopwatch starts.
genotype_matrix <- as.matrix(genotypes[, -"framid"])
storage.mode(genotype_matrix) <- "double"
for (column in seq_len(ncol(genotype_matrix))) {
  if (variants$effect_allele[[column]] == "REFERENCE")
    genotype_matrix[, column] <- 2 - genotype_matrix[, column]
  missing <- !is.finite(genotype_matrix[, column])
  if (any(missing))
    genotype_matrix[missing, column] <- mean(genotype_matrix[, column], na.rm = TRUE)
  genotype_matrix[, column] <- genotype_matrix[, column] * variants$weight[[column]]
}
sets <- lapply(windows$variant_keys, function(value) {
  keys <- strsplit(value, ";", fixed = TRUE)[[1L]]
  match(keys, variants$key)
})

formula_for <- function(model) {
  as.formula(paste0(model,
    " ~ Sex + Age + WBC_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred"))
}

null_for <- function(model) {
  cache <- file.path(options$prepared_dir,
    paste0("r_gmmat_null_", tolower(model), ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  family <- if (model == "BMI") gaussian() else binomial(link = "logit")
  fit <- glmmkin(formula_for(model), data = analysis, kins = relationship,
    id = "framid", family = family, method.optim = "AI", verbose = FALSE)
  saveRDS(fit, cache)
  fit
}

score_state <- function(null, weighted_genotypes) {
  residuals <- null$scaled.residuals
  U <- as.vector(crossprod(weighted_genotypes, residuals))
  if (!is.null(null$P)) {
    V <- crossprod(weighted_genotypes,
      crossprod(null$P, weighted_genotypes))
  } else {
    GSigma_iX <- crossprod(weighted_genotypes, null$Sigma_iX)
    V <- crossprod(weighted_genotypes,
      crossprod(null$Sigma_i, weighted_genotypes)) -
      tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, null$cov))
  }
  list(U = U, V = as.matrix(V))
}

rho <- c(0, 0.25, 0.5, 0.75, 1)
run_test <- function(method, null, columns) {
  state <- score_state(null, genotype_matrix[, columns, drop = FALSE])
  if (method == "burden") {
    score <- sum(state$U)
    variance <- sum(state$V)
    return(list(p = pchisq(score^2 / variance, df = 1,
      lower.tail = FALSE), statistic = score / sqrt(variance)))
  }
  if (method == "skat") {
    return(list(p = GMMAT:::.quad_pval(state$U, state$V,
      method = "davies"), statistic = sum(state$U^2)))
  }
  value <- GMMAT:::.skato_pval(state$U, state$V, rho = rho,
    method = "davies")
  if (method == "skat-o")
    return(list(p = value$p, statistic = value$minp))
  list(p = value$p, statistic = value$minp,
    burden_p = value$Burden.pval, skat_p = value$SKAT.pval)
}

timings <- list()
results <- list()
timing_index <- 0L
result_index <- 0L
for (model in options$models) {
  null <- null_for(model)
  for (method in c("burden", "skat", "skat-o", "suite")) {
    invisible(run_test(method, null, sets[[1L]]))
    for (measurement in seq_len(options$measurements)) {
      gc()
      started <- proc.time()[["elapsed"]]
      value <- NULL
      for (columns in sets) value <- run_test(method, null, columns)
      seconds <- proc.time()[["elapsed"]] - started
      if (!is.finite(value$p)) stop("non-finite benchmark checksum")
      timing_index <- timing_index + 1L
      timings[[timing_index]] <- data.table(runtime = "R-GMMAT", model = model,
        method = method, backend = "R-single-thread", measurement = measurement,
        windows = length(sets), seconds = seconds,
        windows_per_second = length(sets) / seconds)
      cat(sprintf(paste0("R-GMMAT model=%s method=%s measurement=%d ",
        "windows=%d seconds=%.6f windows_per_second=%.3f\n"),
        model, method, measurement, length(sets), seconds,
        length(sets) / seconds))
    }
  }
  for (window_index in seq_along(sets)) {
    value <- run_test("suite", null, sets[[window_index]])
    for (method in c("burden", "skat", "skat-o")) {
      result_index <- result_index + 1L
      p <- if (method == "burden") value$burden_p else
        if (method == "skat") value$skat_p else value$p
      results[[result_index]] <- data.table(runtime = "R-GMMAT", model = model,
        window_id = windows$window_id[[window_index]],
        chromosome = windows$chromosome[[window_index]],
        start = windows$start[[window_index]], end = windows$end[[window_index]],
        genes = windows$genes[[window_index]],
        variants = length(sets[[window_index]]), method = method,
        p_value = p, statistic = value$statistic)
    }
  }
}

fwrite(rbindlist(timings), file.path(options$prepared_dir, "r_timings.csv"))
fwrite(rbindlist(results), file.path(options$prepared_dir, "r_results.csv"))
