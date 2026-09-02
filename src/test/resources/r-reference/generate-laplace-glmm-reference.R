# Deterministic lme4 first-order Laplace reference for GlmmLaplaceTest.
suppressPackageStartupMessages(library(lme4))

groups <- 20L
per_group <- 15L
n <- groups * per_group
row <- 0:(n - 1L)
group <- factor(row %/% per_group)
x <- -1 + 2 * (row %% per_group) / (per_group - 1)
random <- 1.1 * sin(1.3 * (row %/% per_group))
eta <- -0.2 + 0.8 * x + random
probability <- plogis(eta)
state <- 1234567
uniform <- numeric(n)
for (i in seq_len(n)) {
  state <- (48271 * state) %% 2147483647
  uniform[i] <- state / 2147483647
}
y <- as.integer(uniform < probability)

fit <- glmer(y ~ x + (1 | group), family = binomial(), nAGQ = 1,
             control = glmerControl(optimizer = "bobyqa",
                 optCtrl = list(maxfun = 200000)))

cat("R.version=", paste(R.version$major, R.version$minor, sep = "."), "\n", sep = "")
cat("lme4.version=", as.character(packageVersion("lme4")), "\n", sep = "")
cat("intercept=", unname(fixef(fit)[1]), "\n", sep = "")
cat("slope=", unname(fixef(fit)[2]), "\n", sep = "")
cat("randomVariance=", as.numeric(VarCorr(fit)$group[1]), "\n", sep = "")
cat("logLikelihood=", as.numeric(logLik(fit)), "\n", sep = "")
selected <- unname(fitted(fit)[c(1L, 150L, 300L)])
cat("fitted0=", format(selected[1], digits = 17), "\n", sep = "")
cat("fitted149=", format(selected[2], digits = 17), "\n", sep = "")
cat("fitted299=", format(selected[3], digits = 17), "\n", sep = "")
