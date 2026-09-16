# Independent base-R reference for covariance-aware multivariate MR-IVW.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
input <- read.delim(file.path(root, "multivariate-mr-input.tsv"), check.names = FALSE)
correlation <- as.matrix(read.delim(file.path(root,
  "multivariate-mr-correlation.tsv"), header = FALSE))

X <- as.matrix(input[, c("beta_exposure_x1", "beta_exposure_x2")])
Y <- as.matrix(input[, c("beta_outcome_y1", "beta_outcome_y2")])
SY <- as.matrix(input[, c("se_outcome_y1", "se_outcome_y2")])
n <- nrow(X); k <- ncol(X); m <- ncol(Y); p <- k * m
information <- matrix(0, p, p)
rhs <- numeric(p)
for (i in seq_len(n)) {
  covariance <- diag(SY[i, ]) %*% correlation %*% diag(SY[i, ])
  inverse <- solve(covariance)
  design <- kronecker(diag(m), matrix(X[i, ], nrow = 1))
  information <- information + t(design) %*% inverse %*% design
  rhs <- rhs + as.vector(t(design) %*% inverse %*% Y[i, ])
}
coefficient <- as.vector(solve(information, rhs))
coefficient_covariance <- solve(information)
standard_error <- sqrt(diag(coefficient_covariance))
q <- 0
for (i in seq_len(n)) {
  covariance <- diag(SY[i, ]) %*% correlation %*% diag(SY[i, ])
  design <- kronecker(diag(m), matrix(X[i, ], nrow = 1))
  residual <- Y[i, ] - as.vector(design %*% coefficient)
  q <- q + as.numeric(t(residual) %*% solve(covariance, residual))
}
q_df <- m * (n - k)
wald <- function(indices) {
  value <- coefficient[indices]
  covariance <- coefficient_covariance[indices, indices, drop = FALSE]
  as.numeric(t(value) %*% solve(covariance, value))
}
values <- c(
  setNames(coefficient, paste0("beta.", 0:(p - 1))),
  setNames(standard_error, paste0("se.", 0:(p - 1))),
  setNames(as.vector(t(coefficient_covariance)), paste0("covariance.", 0:(p*p - 1))),
  q = q, q_df = q_df, q_p = pchisq(q, q_df, lower.tail = FALSE),
  overall = wald(seq_len(p)), overall_p = pchisq(wald(seq_len(p)), p, lower.tail = FALSE),
  exposure_x1 = wald(c(1, 3)), exposure_x2 = wald(c(2, 4)),
  outcome_y1 = wald(c(1, 2)), outcome_y2 = wald(c(3, 4))
)
options(digits = 17)
writeLines(sprintf("%s=%s", names(values), format(values, scientific = TRUE)),
  file.path(root, "multivariate-mr-reference.properties"))
