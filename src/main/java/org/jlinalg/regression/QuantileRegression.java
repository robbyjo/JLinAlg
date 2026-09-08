/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/**
 * Deterministic quantile regression: fit uses a line-searched smoothed loss;
 * fitExact uses the nonsmoothed pinball linear program. Neither uses sampling.
 */
public final class QuantileRegression {
    private QuantileRegression() { }

    /** Minimizes the nonsmoothed pinball loss by a primal-dual linear program.
     * Predictors explicitly include an intercept column when one is wanted. */
    public static QuantileRegressionResult fitExact(double[] response, double[][] predictors,
                                                    double quantile) {
        return QuantileLinearProgram.solve(response, predictors, quantile).fit();
    }

    /** Exact-loss LP with independent iteration and KKT controls; no smoothing parameter. */
    public static QuantileRegressionResult fitExact(double[] response, double[][] predictors,
            double quantile, QuantileLinearProgram.Options options) {
        return QuantileLinearProgram.solve(response, predictors, quantile, options).fit();
    }

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
        double[] scales = RegressionOptimizer.scaleColumns(x, rows, columns);
        var optimized = RegressionOptimizer.minimize(point -> {
            double[] gradient = gradient(response,x,point,rows,columns,quantile,options.smoothing());
            for(int j=0;j<columns;j++) gradient[j]/=rows;
            return new RegressionOptimizer.Evaluation(
                objective(response,x,point,rows,columns,quantile,options.smoothing())/rows,gradient);
        }, columns,options.maximumIterations(),options.relativeTolerance(),options.initialStep());
        double[] beta = optimized.point();
        double objective=objective(response,x,beta,rows,columns,quantile,options.smoothing());
        double[] fitted = fitted(x, beta, rows, columns), residuals = new double[rows];
        for (int row = 0; row < rows; row++) residuals[row] = response[row] - fitted[row];
        for(int j=0;j<columns;j++) beta[j]/=scales[j];
        return new QuantileRegressionResult(beta, fitted, residuals, objective, quantile, optimized.iterations(), optimized.converged());
    }

    public static QuantileRegressionResult fit(double[] response, double[][] predictors, double quantile) {
        return fit(response, predictors, quantile, QuantileRegressionOptions.defaults());
    }

    private static double[] gradient(double[] y, double[] x, double[] beta, int rows, int columns,
                                     double quantile, double smoothing) {
        double[] result = new double[columns];
        for (int row = 0; row < rows; row++) {
            double residual = y[row] - dot(x, row, beta, columns);
            double scaled = -residual / smoothing;
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
