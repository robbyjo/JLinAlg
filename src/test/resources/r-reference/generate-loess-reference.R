# Regenerates the inline LoessTest reference values with base R stats::loess.
# Usage: Rscript generate-loess-reference.R [output.csv]
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "loess-reference.csv"
x <- c(-2, -1.7, -1.2, -.9, -.4, -.1, .15, .3, .55, .8, 1.05, 1.3,
  1.55, 1.9, 2.2, 2.6, 3, 3.5, 4, 4.6)
y <- sin(x) + .15 * x + c(.02, -.03, .01, .04, -.02, .03, -.01, .02,
  -.04, .01, .03, -.02, .01, -.03, .04, -.01, .02, -.02, .01, -.01)
weights <- c(1, 2, 1, 1, 1.5, 1, 1, 2, 1, 1, 1, 1, 1.5, 1, 1, 1, 2, 1, 1, 1)
control <- loess.control(surface = "direct", statistics = "exact",
  trace.hat = "exact", iterations = 4)
gaussian <- loess(y ~ x, weights = weights, span = .6, degree = 2,
  family = "gaussian", control = control)
contaminated <- y
contaminated[[10]] <- contaminated[[10]] + 3
symmetric <- loess(contaminated ~ x, span = .6, degree = 2,
  family = "symmetric", control = control)
query <- c(-2.2, -1, .2, 1, 2.4, 4.8)
rows <- rbind(
  data.frame(case = "gaussian_fit", index = seq_along(x), x = x,
    value = gaussian$fitted, robust_weight = NA_real_,
    trace_hat = gaussian$trace.hat),
  data.frame(case = "gaussian_predict", index = seq_along(query), x = query,
    value = predict(gaussian, query), robust_weight = NA_real_,
    trace_hat = gaussian$trace.hat),
  data.frame(case = "symmetric_fit", index = seq_along(x), x = x,
    value = symmetric$fitted, robust_weight = symmetric$robust,
    trace_hat = symmetric$trace.hat))
write.csv(rows, output, row.names = FALSE, na = "")
cat(normalizePath(output, mustWork = FALSE), "\n")
