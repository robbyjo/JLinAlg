/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.LeastSquaresSolver;

/** Partially linear semiparametric regression: y = X beta + g(z). */
public final class PartiallyLinearRegression {
    private PartiallyLinearRegression() { }

    public static Result fit(double[] response, double[][] linearPredictors, double[] nonlinearPredictor,
                             double bandwidth, int maximumIterations, BackendPolicy backendPolicy) {
        if (response == null || linearPredictors == null || nonlinearPredictor == null
                || response.length != linearPredictors.length || response.length != nonlinearPredictor.length
                || response.length < 3 || maximumIterations < 1 || !(bandwidth > 0) || !Double.isFinite(bandwidth))
            throw new IllegalArgumentException("partially linear inputs are invalid");
        int rows = response.length, columns = linearPredictors[0] == null ? 0 : linearPredictors[0].length;
        if (columns < 1 || rows <= columns) throw new IllegalArgumentException("linear design is invalid");
        double[] x = flatten(linearPredictors, rows, columns), beta = new double[columns], smooth = new double[rows];
        double[] fitted = new double[rows]; boolean converged = false; int iterations = 0;
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            for (int iteration = 1; iteration <= maximumIterations; iteration++) {
                iterations = iteration;
                double[] adjusted = new double[rows];
                for (int row = 0; row < rows; row++) adjusted[row] = response[row] - smooth[row];
                double[] nextBeta = LeastSquaresSolver.solve(x, adjusted, rows, columns, false, backend).coefficients();
                double[] partial = new double[rows];
                for (int row = 0; row < rows; row++) { double value = 0.0; for (int column = 0; column < columns; column++) value += x[row * columns + column] * nextBeta[column]; partial[row] = response[row] - value; }
                double[] nextSmooth = KernelRegression.predict(nonlinearPredictor, partial, nonlinearPredictor, bandwidth);
                double change = 0.0;
                for (int row = 0; row < rows; row++) change = Math.max(change, Math.abs(nextSmooth[row] - smooth[row]));
                beta = nextBeta; smooth = nextSmooth;
                if (change <= 1e-8 * (1.0 + maxAbs(smooth))) { converged = true; break; }
            }
        }
        for (int row = 0; row < rows; row++) { fitted[row] = smooth[row]; for (int column = 0; column < columns; column++) fitted[row] += x[row * columns + column] * beta[column]; }
        double[] residuals = new double[rows]; for (int row = 0; row < rows; row++) residuals[row] = response[row] - fitted[row];
        return new Result(beta, smooth, fitted, residuals, bandwidth, iterations, converged);
    }

    public static Result fit(double[] response, double[][] linearPredictors, double[] nonlinearPredictor, double bandwidth) {
        return fit(response, linearPredictors, nonlinearPredictor, bandwidth, 100, BackendPolicy.PREFERRED);
    }

    private static double maxAbs(double[] values) { double result = 0.0; for (double value : values) result = Math.max(result, Math.abs(value)); return result; }
    private static double[] flatten(double[][] values, int rows, int columns) {
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) { if (values[row] == null || values[row].length != columns) throw new IllegalArgumentException("linear predictor rows must have equal widths"); for (double value : values[row]) if (!Double.isFinite(value)) throw new IllegalArgumentException("linear predictors must be finite"); System.arraycopy(values[row], 0, result, row * columns, columns); }
        return result;
    }

    /** Semiparametric coefficients, smooth effect, fitted values, and residuals. */
    public record Result(double[] coefficients, double[] smoothEffect, double[] fittedValues,
                         double[] residuals, double bandwidth, int iterations, boolean converged) {
        public Result { coefficients = coefficients.clone(); smoothEffect = smoothEffect.clone(); fittedValues = fittedValues.clone(); residuals = residuals.clone(); }
        public double[] coefficients() { return coefficients.clone(); }
        public double[] beta() { return coefficients(); }
        public double[] smoothEffect() { return smoothEffect.clone(); }
        public double[] fittedValues() { return fittedValues.clone(); }
        public double[] residuals() { return residuals.clone(); }
    }
}
