/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Result of a quantile regression fit. */
public record QuantileRegressionResult(double[] coefficients, double[] fittedValues,
                                       double[] residuals, double objective,
                                       double quantile, int iterations, boolean converged) {
    public QuantileRegressionResult {
        coefficients = coefficients.clone(); fittedValues = fittedValues.clone();
        residuals = residuals.clone();
    }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] beta() { return coefficients(); }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
}
