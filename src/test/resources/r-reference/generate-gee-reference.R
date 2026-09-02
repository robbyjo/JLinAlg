# JLinAlg GEE cross-language reference generator.
# Install geepack 1.3.13 and geer 0.1.0 into JLINALG_R_REFERENCE_LIBRARY.

options(digits = 17, contrasts = c("contr.treatment", "contr.poly"))
arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
script_directory <- dirname(normalizePath(script_argument))
project_root <- normalizePath(file.path(script_directory, "../../../.."))
reference_library <- Sys.getenv(
  "JLINALG_R_REFERENCE_LIBRARY",
  file.path(project_root, ".r-reference-lib")
)
.libPaths(c(reference_library, .libPaths()))

cat("R.version=", R.version.string, "\n", sep = "")
cat("geepack.version=", as.character(packageVersion("geepack")), "\n", sep = "")
cat("geer.version=", as.character(packageVersion("geer")), "\n", sep = "")

emit <- function(name, value) {
  cat(name, "=")
  dput(unname(value))
}

clusters <- 18
id <- repeated <- response <- numeric()
for (cluster in 0:(clusters - 1)) {
  size <- if (cluster < 3) 2 else 3
  cluster_effect <- (cluster %% 5 - 2) * 0.25
  for (visit in 0:(size - 1)) {
    id <- c(id, cluster)
    repeated <- c(repeated, visit)
    response <- c(response, 1.5 + 0.7 * visit + cluster_effect +
      ((cluster + visit) %% 3 - 1) * 0.1)
  }
}
data <- data.frame(response, visit = repeated, id)

geepack_fit <- geepack::geeglm(
  response ~ visit,
  id = id,
  waves = visit,
  data = data,
  family = gaussian(),
  corstr = "exchangeable"
)
emit("gee.gaussian.geepack.coefficients", coef(geepack_fit))
emit("gee.gaussian.geepack.alpha", geepack_fit$geese$alpha)
emit("gee.gaussian.geepack.robust_covariance", as.vector(t(vcov(geepack_fit))))

geer_fit <- geer::geewa(
  response ~ visit,
  id = id,
  repeated = visit,
  data = data,
  family = gaussian(),
  corstr = "exchangeable",
  method = "gee"
)
emit("gee.gaussian.geer.coefficients", coef(geer_fit))
emit("gee.gaussian.geer.alpha", geer_fit$alpha)
emit("gee.gaussian.geer.robust_covariance",
     as.vector(t(vcov(geer_fit, cov_type = "robust"))))
emit("gee.gaussian.geer.bias_corrected_covariance",
     as.vector(t(vcov(geer_fit, cov_type = "bias-corrected"))))

small_sample_data <- subset(data, id >= 3)
kc_fit <- geesmv::GEE.var.kc(
  response ~ visit, id = "id", family = "gaussian",
  data = small_sample_data, corstr = "exchangeable"
)
fg_fit <- geesmv::GEE.var.fg(
  response ~ visit, id = "id", family = "gaussian",
  data = small_sample_data, corstr = "exchangeable"
)
emit("gee.gaussian.kc.variance", kc_fit$cov.beta)
emit("gee.gaussian.fg.variance", fg_fit$cov.beta)

# Mean/scale/correlation estimating equations with a log scale link.
scale_id <- rep(0:29, each = 4)
scale_visit <- rep(0:3, times = 30)
scale_response <- 1.1 + 0.3 * scale_visit +
  exp(0.1 + 0.12 * scale_visit) * sin((seq_along(scale_id) + 2) * 1.7) * 0.18
scale_data <- data.frame(
  response = scale_response,
  visit = scale_visit,
  id = scale_id
)
scale_fit <- geepack::geese(
  response ~ visit,
  sformula = ~ visit,
  id = id,
  waves = visit,
  data = scale_data,
  family = gaussian(),
  corstr = "exchangeable",
  sca.link = "log"
)
emit("gee.scale.coefficients", scale_fit$beta)
emit("gee.scale.gamma", scale_fit$gamma)
emit("gee.scale.alpha", scale_fit$alpha)

# Binomial exchangeable oracle with deterministic pseudo-uniform outcomes.
id <- rep(0:29, each = 4)
visit <- rep(0:3, times = 30)
x <- rep(((0:29) %% 7 - 3) / 3, each = 4)
probability <- plogis(-0.4 + 0.55 * x + 0.25 * visit)
uniform <- (((id * 17 + visit * 13) %% 97) + 0.5) / 97
response <- as.numeric(uniform < probability)
binary_data <- data.frame(response, x, visit, id)
binomial_fit <- geepack::geeglm(
  response ~ x + visit,
  id = id,
  waves = visit,
  data = binary_data,
  family = binomial(),
  corstr = "exchangeable"
)
emit("gee.binomial.coefficients", coef(binomial_fit))
emit("gee.binomial.alpha", binomial_fit$geese$alpha)
emit("gee.binomial.robust_covariance", as.vector(t(vcov(binomial_fit))))

# Poisson AR(1) oracle.
id <- rep(0:23, each = 4)
visit <- rep(0:3, times = 24)
x <- rep(((0:23) %% 6 - 2.5) / 2.5, each = 4)
mean <- exp(0.25 + 0.32 * x + 0.14 * visit)
response <- pmax(0, floor(mean + ((id * 11 + visit * 7) %% 5 - 2) * 0.35))
count_data <- data.frame(response, x, visit, id)
poisson_fit <- geepack::geeglm(
  response ~ x + visit,
  id = id,
  waves = visit,
  data = count_data,
  family = poisson(),
  corstr = "ar1"
)
emit("gee.poisson.coefficients", coef(poisson_fit))
emit("gee.poisson.alpha", poisson_fit$geese$alpha)
emit("gee.poisson.robust_covariance", as.vector(t(vcov(poisson_fit))))

# Weighted, offset, irregular-wave Gaussian oracle.
id <- visit <- response <- weight <- offset <- numeric()
for (cluster in 0:19) {
  visits <- 0:3
  if (cluster %% 4 == 0) visits <- visits[visits != 2]
  for (current in visits) {
    id <- c(id, cluster)
    visit <- c(visit, current)
    weight <- c(weight, 1 + (cluster + current) %% 3 * 0.25)
    offset <- c(offset, 0.08 * current)
    response <- c(response, 0.9 + 0.45 * current + 0.2 * (cluster %% 2) +
      0.08 * current + ((cluster + 2 * current) %% 5 - 2) * 0.04)
  }
}
irregular <- data.frame(response, visit, id, weight, offset)
weighted_fit <- geepack::geeglm(
  response ~ visit,
  id = id,
  waves = visit,
  weights = weight,
  offset = offset,
  data = irregular,
  family = gaussian(),
  corstr = "ar1"
)
emit("gee.weighted_offset.coefficients", coef(weighted_fit))
emit("gee.weighted_offset.alpha", weighted_fit$geese$alpha)
emit("gee.weighted_offset.robust_covariance", as.vector(t(vcov(weighted_fit))))

structure_id <- rep(0:19, each = 4)
structure_visit <- rep(0:3, times = 20)
structure_response <- 1.2 + 0.35 * structure_visit +
  rep(((0:19) %% 5 - 2) * 0.08, each = 4) +
  ((structure_id + structure_visit) %% 3 - 1) * 0.03
structure_data <- data.frame(
  response = structure_response,
  visit = structure_visit,
  id = structure_id
)
for (structure in c("independence", "exchangeable", "ar1", "unstructured")) {
  structure_fit <- geepack::geeglm(
    response ~ visit,
    id = id,
    waves = visit,
    data = structure_data,
    family = gaussian(),
    corstr = structure
  )
  emit(paste0("gee.structure.", structure, ".coefficients"),
       coef(structure_fit))
}
