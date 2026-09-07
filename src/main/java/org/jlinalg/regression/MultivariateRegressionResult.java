/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Multivariate OLS coefficients and residual covariance. */
public record MultivariateRegressionResult(
        double[] coefficients, double[] residualCovariance,
        double[] fittedValues, double[] residuals,
        int observations, int predictorCount, int outcomeCount) {
    public MultivariateRegressionResult {
        coefficients = coefficients.clone(); residualCovariance = residualCovariance.clone();
        fittedValues = fittedValues.clone(); residuals = residuals.clone();
    }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] beta() { return coefficients(); }
    public double[] residualCovariance() { return residualCovariance.clone(); }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
}
