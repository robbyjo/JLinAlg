# Independent ACAT-V/ACAT-O fixtures following Liu et al. (2019) and the
# yaowuliu/ACAT reference implementation. CompQuadForm supplies Davies' tail
# for the two SKAT components; all other calculations use base R.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1) args[[1]] else
  "src/test/resources/acat/reference.tsv"
if (length(args) >= 2) .libPaths(c(args[[2]], .libPaths()))
suppressPackageStartupMessages(library(CompQuadForm))

acat <- function(p, weights = rep(1, length(p))) {
  stopifnot(length(p) > 0, length(p) == length(weights),
            all(is.finite(p)), all(p >= 0), all(p <= 1),
            all(is.finite(weights)), all(weights >= 0), sum(weights) > 0)
  if (any(p == 0) && any(p == 1)) stop("opposing endpoint p-values")
  if (any(p == 0)) return(0)
  if (any(p == 1)) return(1)
  weights <- weights / sum(weights)
  small <- p < 1e-15
  statistic <- sum(weights[small] / (pi * p[small])) +
    sum(weights[!small] * tan((0.5 - p[!small]) * pi))
  if (statistic > 1e15) return(1 / (pi * statistic))
  pcauchy(statistic, lower.tail = FALSE)
}

burden_p <- function(u, v, weights, selected = rep(TRUE, length(u))) {
  weighted_u <- weights[selected] * u[selected]
  weighted_v <- outer(weights[selected], weights[selected]) *
    v[selected, selected, drop = FALSE]
  2 * pnorm(-abs(sum(weighted_u) / sqrt(sum(weighted_v))))
}

skat_p <- function(u, v, weights) {
  q <- sum((weights * u)^2)
  lambda <- eigen(outer(weights, weights) * v,
                  symmetric = TRUE, only.values = TRUE)$values
  lambda <- lambda[lambda > max(abs(lambda)) * 1e-12]
  CompQuadForm::davies(q, lambda)$Qq
}

acat_v <- function(u, v, maf, mac, shapes, threshold = 10) {
  coefficients <- dbeta(maf, shapes[[1]], shapes[[2]])
  p_weights <- (coefficients / dbeta(maf, 0.5, 0.5))^2
  collapsed <- mac <= threshold
  p <- numeric(0)
  w <- numeric(0)
  if (any(collapsed)) {
    p <- c(p, burden_p(u, v, coefficients, collapsed))
    mean_maf <- mean(maf[collapsed])
    w <- c(w, (dbeta(mean_maf, shapes[[1]], shapes[[2]]) /
               dbeta(mean_maf, 0.5, 0.5))^2)
  }
  dense <- which(!collapsed)
  if (length(dense)) {
    p <- c(p, 2 * pnorm(-abs(u[dense] / sqrt(diag(v)[dense]))))
    w <- c(w, p_weights[dense])
  }
  acat(p, w)
}

evaluate <- function(id, u, v, maf, mac) {
  w_1_25 <- dbeta(maf, 1, 25)
  w_1_1 <- dbeta(maf, 1, 1)
  values <- c(
    skat_p(u, v, w_1_25),
    skat_p(u, v, w_1_1),
    burden_p(u, v, w_1_25),
    burden_p(u, v, w_1_1),
    acat_v(u, v, maf, mac, c(1, 25)),
    acat_v(u, v, maf, mac, c(1, 1)))
  data.frame(
    id = id,
    u = paste(format(u, digits = 17, scientific = FALSE), collapse = ","),
    v = paste(format(c(t(v)), digits = 17, scientific = FALSE), collapse = ","),
    maf = paste(format(maf, digits = 17, scientific = FALSE), collapse = ","),
    mac = paste(format(mac, digits = 17, scientific = FALSE), collapse = ","),
    skat_1_25 = values[[1]], skat_1_1 = values[[2]],
    burden_1_25 = values[[3]], burden_1_1 = values[[4]],
    acat_v_1_25 = values[[5]], acat_v_1_1 = values[[6]],
    acat_o = acat(values), check.names = FALSE)
}

fixtures <- list(
  evaluate("mixed", c(1.8, -0.7, 2.1, 0.3),
    matrix(c(4, .4, .2, 0, .4, 3, .3, .1,
             .2, .3, 2, .2, 0, .1, .2, 1.5), 4, 4, byrow = TRUE),
    c(.0005, .002, .01, .03), c(2, 8, 40, 120)),
  evaluate("dense", c(-1.1, 2.4, .8),
    matrix(c(2.5, .2, -.1, .2, 3.2, .4, -.1, .4, 1.8),
           3, 3, byrow = TRUE),
    c(.003, .0075, .02), c(12, 30, 80)),
  evaluate("collapsed", c(.4, -1.2, 1.7),
    matrix(c(1.2, .1, .05, .1, 1.8, .2, .05, .2, 2.1),
           3, 3, byrow = TRUE),
    c(.00025, .001, .00225), c(1, 4, 9)),
  evaluate("correlated", c(3.7, -2.2, 1.1, 2.9, -.6),
    crossprod(matrix(c(1, .2, 0, 0, 0,
                       .1, 1, .3, 0, 0,
                       0, .2, 1, .2, 0,
                       0, 0, .1, 1, .25,
                       0, 0, 0, .15, 1), 5, 5, byrow = TRUE)),
    c(.0004, .0014, .004, .012, .04), c(2, 7, 20, 60, 200)))

result <- do.call(rbind, fixtures)
options(digits = 17)
write.table(result, output, sep = "\t", row.names = FALSE, quote = FALSE,
            na = "NA")
cat(sprintf("generic_acat=%.17g\n", acat(c(.02, .0004, .2, .1, .8))))
