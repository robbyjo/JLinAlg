/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Gaussian multivariate ordinary least squares with shared predictors. */
public final class MultivariateRegression {
    private MultivariateRegression() { }

    public static MultivariateRegressionResult fit(double[][] response,
            double[][] predictors, BackendPolicy backendPolicy) {
        if (response == null || predictors == null || response.length < 2
                || predictors.length != response.length || response[0] == null
                || predictors[0] == null)
            throw new IllegalArgumentException("response and predictor rows are required");
        int rows = response.length, outcomes = response[0].length;
        int columns = predictors[0].length;
        if (outcomes < 1 || columns < 1 || rows <= columns)
            throw new IllegalArgumentException("multivariate regression needs full-rank rows");
        for(double[] row:response) if(row==null || row.length!=outcomes)
            throw new IllegalArgumentException("response rows must have equal widths");
        double[] x = flatten(predictors, rows, columns);
        double[] fitted = new double[rows * outcomes];
        double[] coefficients = new double[columns * outcomes];
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            var qr=backend.dgeqp3(x,rows,columns);
            if(qr.rank()!=columns)throw new IllegalArgumentException("multivariate design is rank deficient");
            for (int outcome = 0; outcome < outcomes; outcome++) {
                double[] y = column(response, outcome, rows);
                double[] solved = qr.solveLeastSquares(y);
                System.arraycopy(solved, 0, coefficients,
                    outcome * columns, columns);
                for (int row = 0; row < rows; row++) {
                    double value = 0.0;
                    for (int column = 0; column < columns; column++)
                        value += x[row * columns + column]
                            * solved[column];
                    fitted[row * outcomes + outcome] = value;
                }
            }
        }
        double[] residualCovariance = new double[outcomes * outcomes];
        for (int row = 0; row < rows; row++) for (int left = 0; left < outcomes; left++) {
            double leftResidual = response[row][left] - fitted[row * outcomes + left];
            for (int right = 0; right < outcomes; right++)
                residualCovariance[left * outcomes + right] += leftResidual
                    * (response[row][right] - fitted[row * outcomes + right]);
        }
        double divisor = rows - columns;
        for (int i = 0; i < residualCovariance.length; i++) residualCovariance[i] /= divisor;
        return new MultivariateRegressionResult(coefficients, residualCovariance,
            fitted, residuals(response, fitted), rows, columns, outcomes);
    }

    public static MultivariateRegressionResult fit(double[][] response,
            double[][] predictors) {
        return fit(response, predictors, BackendPolicy.PREFERRED);
    }

    private static double[] flatten(double[][] values, int rows, int columns) {
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            if (values[row] == null || values[row].length != columns)
                throw new IllegalArgumentException("rows must have equal widths");
            for (double value : values[row]) if (!Double.isFinite(value))
                throw new IllegalArgumentException("predictors must be finite");
            System.arraycopy(values[row], 0, result, row * columns, columns);
        }
        return result;
    }
    private static double[] column(double[][] values, int column, int rows) {
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) {
            if (values[row] == null || values[row].length <= column
                    || !Double.isFinite(values[row][column]))
                throw new IllegalArgumentException("response rows must be finite and equal width");
            result[row] = values[row][column];
        }
        return result;
    }
    private static double[] residuals(double[][] response, double[] fitted) {
        int rows = response.length, outcomes = response[0].length;
        double[] result = new double[fitted.length];
        for (int row = 0; row < rows; row++) for (int outcome = 0; outcome < outcomes; outcome++)
            result[row * outcomes + outcome] = response[row][outcome]
                - fitted[row * outcomes + outcome];
        return result;
    }
}
