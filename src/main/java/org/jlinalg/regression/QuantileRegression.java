/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/**
 * Quantile regression using a deterministic, line-searched smoothed pinball loss.
 * This is a stable first implementation and does not use bootstrap sampling.
 */
public final class QuantileRegression {
    private QuantileRegression() { }

    public static QuantileRegressionResult fit(double[] response, double[][] predictors,
                                               double quantile,
                                               QuantileRegressionOptions options) {
        if (response == null || predictors == null || options == null
                || response.length != predictors.length || response.length < 2
                || !(quantile > 0) || !(quantile < 1) || !Double.isFinite(quantile))
            throw new IllegalArgumentException("quantile inputs are invalid");
        int rows = response.length, columns = predictors[0] == null ? 0 : predictors[0].length;
        if (columns < 1 || rows <= columns) throw new IllegalArgumentException("quantile design is invalid");
        double[] x = flatten(predictors, rows, columns);
        for (double value : response) if (!Double.isFinite(value)) throw new IllegalArgumentException("response must be finite");
        double[] beta = new double[columns];
        double objective = objective(response, x, beta, rows, columns, quantile, options.smoothing());
        boolean converged = false; int iterations = 0;
        for (int iteration = 1; iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            double[] gradient = gradient(response, x, beta, rows, columns, quantile, options.smoothing());
            double step = options.initialStep();
            double[] candidate = beta.clone();
            double candidateObjective = Double.POSITIVE_INFINITY;
            for (int attempt = 0; attempt < 50; attempt++) {
                for (int column = 0; column < columns; column++) candidate[column] = beta[column] - step * gradient[column] / rows;
                candidateObjective = objective(response, x, candidate, rows, columns, quantile, options.smoothing());
                if (candidateObjective <= objective) break;
                step *= 0.5;
            }
            double change = 0.0;
            for (int column = 0; column < columns; column++) change = Math.max(change,
                Math.abs(candidate[column] - beta[column]) / (1.0 + Math.abs(beta[column])));
            beta = candidate; double improvement = objective - candidateObjective; objective = candidateObjective;
            if (change <= options.relativeTolerance() && Math.abs(improvement) <= options.relativeTolerance()) { converged = true; break; }
        }
        double[] fitted = fitted(x, beta, rows, columns), residuals = new double[rows];
        for (int row = 0; row < rows; row++) residuals[row] = response[row] - fitted[row];
        return new QuantileRegressionResult(beta, fitted, residuals, objective, quantile, iterations, converged);
    }

    public static QuantileRegressionResult fit(double[] response, double[][] predictors, double quantile) {
        return fit(response, predictors, quantile, QuantileRegressionOptions.defaults());
    }

    private static double[] gradient(double[] y, double[] x, double[] beta, int rows, int columns,
                                     double quantile, double smoothing) {
        double[] result = new double[columns];
        for (int row = 0; row < rows; row++) {
            double residual = y[row] - dot(x, row, beta, columns);
            double scaled = residual / smoothing;
            double probability = scaled >= 0.0
                ? 1.0 / (1.0 + Math.exp(-scaled))
                : Math.exp(scaled) / (1.0 + Math.exp(scaled));
            double derivative = -(quantile - probability);
            for (int column = 0; column < columns; column++) result[column] += derivative * x[row * columns + column];
        }
        return result;
    }
    private static double objective(double[] y, double[] x, double[] beta, int rows, int columns,
                                    double quantile, double smoothing) {
        double result = 0.0;
        for (int row = 0; row < rows; row++) {
            double residual = y[row] - dot(x, row, beta, columns);
            result += quantile * residual + smoothing * log1pExp(-residual / smoothing);
        }
        return result;
    }
    private static double dot(double[] x, int row, double[] beta, int columns) {
        double result = 0.0; for (int column = 0; column < columns; column++) result += x[row * columns + column] * beta[column]; return result;
    }
    private static double log1pExp(double value) {
        return value > 0.0 ? value + Math.log1p(Math.exp(-value))
            : Math.log1p(Math.exp(value));
    }
    private static double[] fitted(double[] x, double[] beta, int rows, int columns) {
        double[] result = new double[rows]; for (int row = 0; row < rows; row++) result[row] = dot(x, row, beta, columns); return result;
    }
    private static double[] flatten(double[][] values, int rows, int columns) {
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            if (values[row] == null || values[row].length != columns) throw new IllegalArgumentException("predictor rows must have equal widths");
            for (double value : values[row]) if (!Double.isFinite(value)) throw new IllegalArgumentException("predictors must be finite");
            System.arraycopy(values[row], 0, result, row * columns, columns);
        }
        return result;
    }
}
