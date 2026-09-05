# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later
suppressPackageStartupMessages(library(glmmTMB))

set.seed(20260905)
groups <- 18L
per_group <- 25L
n <- groups * per_group
group <- factor(rep(seq_len(groups), each = per_group))
x <- rep(seq(-1, 1, length.out = per_group), groups)
bc <- rnorm(groups, sd = 0.55)
bz <- rnorm(groups, sd = 0.7)
eta_count <- 0.55 + 0.35 * x + bc[group]
eta_zero <- -0.9 - 0.25 * x + bz[group]
pi <- plogis(eta_zero)
mu <- exp(eta_count)
structural <- runif(n) < pi
y_zip <- ifelse(structural, 0L, rpois(n, mu))
y_zinb <- ifelse(structural, 0L, rnbinom(n, mu = mu, size = 2.4))
data <- data.frame(y_zip, y_zinb, x, group)

zip <- glmmTMB(y_zip ~ x + (1 | group), ziformula = ~ x + (1 | group),
  family = poisson, data = data)
zinb <- glmmTMB(y_zinb ~ x + (1 | group),
  ziformula = ~ x, family = nbinom2, data = data)

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else
  file.path("src", "test", "resources", "r-reference")
write.table(data, file.path(root, "zero-inflated-mixed-data.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t")

values <- c(
  generator = "glmmTMB 1.1.14 under R 4.6.1",
  zip_count_intercept = fixef(zip)$cond[[1L]],
  zip_count_x = fixef(zip)$cond[[2L]],
  zip_zero_intercept = fixef(zip)$zi[[1L]],
  zip_zero_x = fixef(zip)$zi[[2L]],
  zip_count_sd = attr(VarCorr(zip)$cond$group, "stddev")[[1L]],
  zip_zero_sd = attr(VarCorr(zip)$zi$group, "stddev")[[1L]],
  zip_count_intercept_se = summary(zip)$coefficients$cond[1L, 2L],
  zip_count_x_se = summary(zip)$coefficients$cond[2L, 2L],
  zip_zero_intercept_se = summary(zip)$coefficients$zi[1L, 2L],
  zip_zero_x_se = summary(zip)$coefficients$zi[2L, 2L],
  zip_log_likelihood = as.numeric(logLik(zip)),
  zinb_count_intercept = fixef(zinb)$cond[[1L]],
  zinb_count_x = fixef(zinb)$cond[[2L]],
  zinb_zero_intercept = fixef(zinb)$zi[[1L]],
  zinb_zero_x = fixef(zinb)$zi[[2L]],
  zinb_size = sigma(zinb),
  zinb_count_sd = attr(VarCorr(zinb)$cond$group, "stddev")[[1L]],
  zinb_log_likelihood = as.numeric(logLik(zinb))
)
connection <- file(file.path(root, "zero-inflated-mixed-glmmtmb.properties"),
  open = "wt")
on.exit(close(connection))
for (name in names(values)) writeLines(paste0(name, "=", values[[name]]), connection)
