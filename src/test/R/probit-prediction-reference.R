# Regenerates ProbitPredictionRReferenceTest with base R only.
options(digits = 17)
x <- c(-2, -1.7, -1.4, -1.1, -.8, -.5, -.2, .1,
       .4, .7, 1, 1.3, 1.6, 1.9, 2.2, 2.5)
z <- rep(c(0, 1), 8)
y <- c(0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1)
X <- cbind(1, x, z)
fit <- glm(y ~ x + z, family = binomial(link = "probit"),
           control = glm.control(epsilon = 1e-12, maxit = 100))
print(coef(fit))
print(sqrt(diag(vcov(fit))))
print(c(deviance = deviance(fit), logLik = as.numeric(logLik(fit)),
        aic = AIC(fit)))

new_design <- rbind(c(1, -1, 0), c(1, 0, 1), c(1, 1, 0))
prediction <- predict(fit,
    newdata = data.frame(x = new_design[, 2], z = new_design[, 3]),
    type = "link", se.fit = TRUE)
eta <- prediction$fit
link_se <- prediction$se.fit
print(cbind(eta, link_se, mean = pnorm(eta),
    mean_se = dnorm(eta) * link_se,
    lower = pnorm(eta - qnorm(.975) * link_se),
    upper = pnorm(eta + qnorm(.975) * link_se)))

V <- vcov(fit)
beta <- coef(fit)
scenario <- function(design) {
    eta <- drop(design %*% beta)
    list(mean = mean(pnorm(eta)), gradient = colMeans(dnorm(eta) * design))
}
first <- scenario(cbind(1, x, 1))
second <- scenario(cbind(1, x, 0))
covariance <- drop(t(first$gradient) %*% V %*% second$gradient)
difference_gradient <- first$gradient - second$gradient
difference <- first$mean - second$mean
difference_se <- sqrt(drop(t(difference_gradient) %*% V %*%
    difference_gradient))
log_ratio_gradient <- first$gradient / first$mean -
    second$gradient / second$mean
log_ratio_se <- sqrt(drop(t(log_ratio_gradient) %*% V %*%
    log_ratio_gradient))
ratio <- first$mean / second$mean
print(c(first = first$mean, second = second$mean, covariance = covariance,
    difference = difference, difference_se = difference_se,
    ratio = ratio, ratio_se = ratio * log_ratio_se,
    log_ratio_se = log_ratio_se))

column <- 2
eta <- drop(X %*% beta)
effect <- mean(dnorm(eta) * beta[column])
gradient <- colMeans((-eta * dnorm(eta) * beta[column]) * X)
gradient[column] <- gradient[column] + mean(dnorm(eta))
effect_se <- sqrt(drop(t(gradient) %*% V %*% gradient))
print(c(average_marginal_effect = effect, standard_error = effect_se))
