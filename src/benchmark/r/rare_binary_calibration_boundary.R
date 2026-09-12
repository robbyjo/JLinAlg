# Independent finite-support example explaining why normal U/V tails are not
# binary rare-case calibration. Base R only; no sampling or external packages.
# Run from the repository root. This is an assessment fixture, not an estimator.
n <- 1000L
cases <- 1L
carriers <- 10L
observed <- 1L
p <- cases / n
u <- observed - carriers * p
# Efficient logistic information after projecting out the intercept.
v <- p * (1-p) * (carriers - carriers^2/n)
support <- seq.int(max(0L, cases-(n-carriers)), min(cases, carriers))
mass <- dhyper(support, cases, n-cases, carriers)
exact_abs_score_tail <- sum(mass[abs(support-carriers*p) >= abs(u)-1e-14])
stopifnot(abs(exact_abs_score_tail - .01) < 1e-14)
normal_tail <- 2 * pnorm(-abs(u)/sqrt(v))
stopifnot(normal_tail < 1e-20)
result <- data.frame(n=n, cases=cases, carriers=carriers, observed=observed,
                     score=u, logistic_information=v,
                     normal_p=normal_tail, conditional_abs_score_p=exact_abs_score_tail)
options(digits=17)
write.table(result, "src/test/resources/raremetal/binary-calibration-boundary.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
print(result)
