/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import org.jlinalg.internal.MatrixOps;

/** Immutable Gaussian ridge/LASSO/elastic-net fit at one penalty value. */
public final class PenalizedRegressionResult {
    private final double intercept;
    private final double[] coefficients;
    private volatile double[] fittedValues;
    private volatile double[] residuals;
    private final double[] response;
    private final double[] predictors;
    private final int rows;
    private final int columns;
    private final double lambda;
    private final double alpha;
    private final double objective;
    private final double weightedResidualSumOfSquares;
    private final int activeCoefficientCount;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;

    PenalizedRegressionResult(
            double intercept,
            double[] coefficients,
            double[] response,
            double[] predictors,
            int rows,
            int columns,
            double lambda,
            double alpha,
            double objective,
            double weightedResidualSumOfSquares,
            int activeCoefficientCount,
            int iterations,
            boolean converged,
            String convergenceMessage) {
        this.intercept = intercept;
        this.coefficients = coefficients.clone();
        this.response = response;
        this.predictors = predictors;
        this.rows = rows;
        this.columns = columns;
        this.lambda = lambda;
        this.alpha = alpha;
        this.objective = objective;
        this.weightedResidualSumOfSquares = weightedResidualSumOfSquares;
        this.activeCoefficientCount = activeCoefficientCount;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
    }

    public double intercept() { return intercept; }
    public double[] coefficients() { return coefficients.clone(); }
    /** Alias for coefficients, consistent with the association-oriented APIs. */
    public double[] beta() { return coefficients(); }
    public double[] fittedValues() {
        materializePredictions();
        return fittedValues.clone();
    }
    public double[] residuals() {
        materializePredictions();
        return residuals.clone();
    }
    public double lambda() { return lambda; }
    public double alpha() { return alpha; }
    public double objective() { return objective; }
    public double weightedResidualSumOfSquares() {
        return weightedResidualSumOfSquares;
    }
    public int activeCoefficientCount() { return activeCoefficientCount; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }

    private void materializePredictions() {
        if (fittedValues != null) return;
        synchronized (this) {
            if (fittedValues != null) return;
            double[] fitted = new double[rows];
            double[] errors = new double[rows];
            for (int row = 0; row < rows; row++) {
                double value = intercept;
                for (int column = 0; column < columns; column++) {
                    value += predictors[row * columns + column]
                        * coefficients[column];
                }
                fitted[row] = value;
                errors[row] = response[row] - value;
            }
            residuals = errors;
            fittedValues = fitted;
        }
    }

    /** Predicts on a conventional matrix in the original predictor scale. */
    public double[] predict(double[][] predictors) {
        if (predictors == null) {
            throw new IllegalArgumentException("predictors are required");
        }
        double[] values = MatrixOps.rowMajor(predictors, predictors.length);
        int columns = predictors[0].length;
        if (columns != coefficients.length) {
            throw new IllegalArgumentException(
                "prediction columns must equal the fitted coefficient count");
        }
        double[] result = new double[predictors.length];
        for (int row = 0; row < predictors.length; row++) {
            double value = intercept;
            for (int column = 0; column < columns; column++) {
                value += values[row * columns + column] * coefficients[column];
            }
            result[row] = value;
        }
        return result;
    }
}
