# SPDX-License-Identifier: GPL-2.0-or-later
# Independent base-R reference for individual-level 2SLS. No Java output is
# consumed. The calculations use structural residuals and the GMM/IV sandwich,
# not lm(y ~ fitted.endogenous) standard errors.

n <- 180L
index <- 0:(n - 1L)
cluster <- index %/% 6L
w <- ((index %% 13L) - 6) / 4 + 0.2 * sin((index + 1) * 0.17)
z1 <- sin((index + 1) * 0.31) + 0.10 * w
z2 <- cos((index + 2) * 0.19) - 0.15 * w
z3 <- ((index * 7L) %% 17L - 8) / 6 + 0.2 * sin((index + 3) * 0.07)
cluster_effect <- ((cluster %% 7L) - 3) * 0.12
base_error <- sin((index + 1) * 1.11) + 0.4 * cos((index + 2) * 0.47)
u <- (0.35 + 0.04 * (index %% 5L)) * base_error + cluster_effect
v1 <- 0.25 * cos((index + 1) * 0.73) + 0.10 * sin((index + 2) * 0.13)
v2 <- 0.22 * sin((index + 1) * 0.61) - 0.08 * cos((index + 2) * 0.29)
x1 <- 0.70 * z1 + 0.35 * z2 + 0.25 * w + 0.65 * u + v1
x2 <- 0.55 * z2 - 0.40 * z3 - 0.20 * w - 0.45 * u + v2
y <- 1.10 + 0.40 * w + 1.25 * x1 - 0.75 * x2 + u

W <- cbind(intercept = 1, w = w)
X <- cbind(x1 = x1, x2 = x2)
Z <- cbind(z1 = z1, z2 = z2, z3 = z3)
A <- cbind(W, X)
B <- cbind(W, Z)
p <- ncol(A)
l <- ncol(B)
q <- ncol(Z)
groups <- length(unique(cluster))

# Pivoted QR projection; fitted() avoids constructing the n by n projector.
first_stage <- qr(B, LAPACK = TRUE)
Pi <- qr.coef(first_stage, X)
Xhat <- B %*% Pi
Ahat <- cbind(W, Xhat)
beta <- drop(solve(crossprod(Ahat), crossprod(Ahat, y)))
bread <- solve(crossprod(Ahat))
residual <- drop(y - A %*% beta)
rss <- sum(residual^2)

outer_meat <- function(design, errors) {
  scores <- design * errors
  crossprod(scores)
}

cluster_meat <- function(design, errors) {
  scores <- rowsum(design * errors, cluster, reorder = FALSE)
  crossprod(scores)
}

sandwich <- function(bread, meat) bread %*% meat %*% bread

covariances <- list(
  homoskedastic = bread * rss / (n - p),
  hc0 = sandwich(bread, outer_meat(Ahat, residual)),
  hc1 = sandwich(bread, outer_meat(Ahat, residual)) * n / (n - p),
  cluster_cr0 = sandwich(bread, cluster_meat(Ahat, residual)),
  cluster_cr1 = sandwich(bread, cluster_meat(Ahat, residual)) *
    groups / (groups - 1) * (n - 1) / (n - p)
)

first_stage_diagnostics <- function(column, covariance_name) {
  target <- X[, column]
  fitted_unrestricted <- Xhat[, column]
  errors <- target - fitted_unrestricted
  rss_u <- sum(errors^2)
  rss_r <- sum((target - W %*% qr.coef(qr(W, LAPACK = TRUE), target))^2)
  improvement <- max(0, rss_r - rss_u)
  partial_r2 <- improvement / rss_r
  classical_f <- (improvement / q) / (rss_u / (n - l))
  coefficients <- Pi[, column]
  first_bread <- solve(crossprod(B))
  if (covariance_name == "homoskedastic") {
    wald <- q * classical_f
    effective_f <- classical_f
    p_value <- pf(classical_f, q, n - l, lower.tail = FALSE)
  } else {
    if (covariance_name == "hc0") {
      covariance <- sandwich(first_bread, outer_meat(B, errors))
    } else if (covariance_name == "hc1") {
      covariance <- sandwich(first_bread, outer_meat(B, errors)) * n / (n - l)
    } else if (covariance_name == "cluster_cr0") {
      covariance <- sandwich(first_bread, cluster_meat(B, errors))
    } else {
      covariance <- sandwich(first_bread, cluster_meat(B, errors)) *
        groups / (groups - 1) * (n - 1) / (n - l)
    }
    selected <- (ncol(W) + 1):l
    wald <- drop(t(coefficients[selected]) %*%
      solve(covariance[selected, selected], coefficients[selected]))
    effective_f <- wald / q
    p_value <- if (startsWith(covariance_name, "cluster")) {
      pf(effective_f, q, groups - l, lower.tail = FALSE)
    } else {
      pchisq(wald, q, lower.tail = FALSE)
    }
  }
  c(partial_r2 = partial_r2, classical_f = classical_f,
    wald = wald, effective_f = effective_f, p_value = p_value)
}

rows <- list()
add <- function(metric, values) {
  rows[[length(rows) + 1L]] <<- data.frame(
    metric = metric, index = seq_along(values) - 1L,
    value = as.numeric(values))
}
add("beta", beta)
add("rss", rss)
add("residual_variance", rss / (n - p))
add("instrumented_selected", as.vector(t(Xhat[c(1, 2, 18, 79, 180), ])))
Ahat_scaled <- sweep(Ahat, 2, sqrt(colSums(A^2)), "/")
singular <- svd(Ahat_scaled, nu = 0, nv = 0)$d
add("projected_condition", max(singular) / min(singular))
for (name in names(covariances)) {
  add(paste0(name, "_cov"), as.vector(t(covariances[[name]])))
  for (column in seq_len(ncol(X))) {
    add(paste0(name, "_first_stage_", column - 1L),
      first_stage_diagnostics(column, name))
  }
}
reference <- do.call(rbind, rows)

root <- file.path("src", "test", "resources", "r-reference")
write.csv(data.frame(cluster, w, z1, z2, z3, x1, x2, y),
  file.path(root, "instrumental-variables-data.csv"),
  row.names = FALSE, quote = FALSE)
write.csv(reference,
  file.path(root, "instrumental-variables-reference.csv"),
  row.names = FALSE, quote = FALSE)
writeLines(c(
  paste0("runtime=", R.version.string),
  "implementation=base R pivoted QR and explicit IV/GMM sandwiches",
  "packages=base,stats"
), file.path(root, "instrumental-variables-reference.properties"))

cat(R.version.string, "\n")
cat("wrote", nrow(reference), "reference values\n")
