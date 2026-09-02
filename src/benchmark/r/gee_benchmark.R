# Comparable geepack benchmark for org.jlinalg.benchmark.GeeBenchmark.
arguments <- commandArgs(trailingOnly = TRUE)
clusters <- if (length(arguments) >= 1) as.integer(arguments[[1]]) else 10000L
size <- if (length(arguments) >= 2) as.integer(arguments[[2]]) else 5L

id <- rep(0:(clusters - 1L), each = size)
visit <- rep(0:(size - 1L), times = clusters)
treatment <- id %% 2L
row <- seq_along(id)
cluster_effect <- ((id %% 17L) - 8L) * 0.025
noise <- sin(row * 12.9898) * 0.2
response <- 1.0 + 0.4 * visit + 0.25 * treatment + cluster_effect + noise
data <- data.frame(response, visit, treatment, id)

cat(sprintf("geepack=%s rows=%d clusters=%d\n",
            as.character(packageVersion("geepack")), nrow(data), clusters))
for (structure in c("independence", "exchangeable", "ar1", "unstructured")) {
  elapsed <- system.time({
    fit <- geepack::geeglm(
      response ~ visit + treatment,
      id = id,
      waves = visit,
      data = data,
      family = gaussian(),
      corstr = structure
    )
  })[["elapsed"]]
  cat(sprintf("%s seconds=%.6f beta=%s alpha=%s\n",
              structure, elapsed,
              paste(format(coef(fit), digits = 12), collapse = ","),
              paste(format(fit$geese$alpha, digits = 12), collapse = ",")))
}
