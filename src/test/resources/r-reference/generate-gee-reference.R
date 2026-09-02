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
