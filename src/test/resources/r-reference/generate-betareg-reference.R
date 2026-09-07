# Generates the reference values asserted by BetaRegressionRReferenceTest.
# Requires the CRAN betareg package. The data files are package examples.
library(betareg)

data("GasolineYield", package = "betareg")
data("FoodExpenditure", package = "betareg")

fits <- list(
  gasoline = betareg(yield ~ batch + temp, data = GasolineYield),
  gasoline_variable = betareg(
    yield ~ batch + temp | temp, data = GasolineYield),
  food = betareg(
    I(food / income) ~ income + persons, data = FoodExpenditure)
)

for (name in names(fits)) {
  fit <- fits[[name]]
  cat(name, "\n")
  print(coef(fit), digits = 16)
  print(sqrt(diag(vcov(fit))), digits = 16)
  print(logLik(fit), digits = 16)
}
