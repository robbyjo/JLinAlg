# Generates the deterministic glmmTMB accuracy fixture used by
# BetaMixedModelRReferenceTest. Run from the repository root.
.libPaths(c("build/r-library", .libPaths()))
suppressPackageStartupMessages(library(glmmTMB))

groups <- 40L
per_group <- 15L
n <- groups * per_group
index <- seq_len(n)
group_index <- rep(seq_len(groups), each = per_group)
group <- factor(sprintf("g%02d", group_index))
x <- (((index * 37L) %% 101L) - 50) / 50
u <- 0.70 * sin(seq_len(groups) * 1.7) +
  0.20 * cos(seq_len(groups) * 0.43)
mu <- plogis(-0.35 + 0.90 * x + u[group_index])
generating_phi <- 18
probability <- (((index * 7919L) %% 997L) + 0.5) / 997
y <- qbeta(probability, mu * generating_phi, (1 - mu) * generating_phi)
data <- data.frame(y = y, x = x, group = group)

fit <- glmmTMB(y ~ x + (1 | group), data = data,
  family = beta_family(link = "logit"))

write.table(data, "src/test/resources/r-reference/beta-mixed-glmmtmb.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)
reference <- c(
  package_version = as.character(packageVersion("glmmTMB")),
  beta_0 = unname(fixef(fit)$cond[1]),
  beta_1 = unname(fixef(fit)$cond[2]),
  precision = exp(unname(fixef(fit)$disp[1])),
  variance = unname(VarCorr(fit)$cond$group[1]),
  log_likelihood = as.numeric(logLik(fit))
)
writeLines(sprintf("%s=%s", names(reference),
  format(reference, digits = 17, trim = TRUE)),
  "src/test/resources/r-reference/beta-mixed-glmmtmb.properties")
print(reference)
