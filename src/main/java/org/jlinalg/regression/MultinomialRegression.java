/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Baseline-category multinomial logistic regression by monotone gradient ascent. */
public final class MultinomialRegression {
    private MultinomialRegression() { }

    public static MultinomialRegressionResult fit(int[] response,
            double[][] predictors, int classCount, MultinomialOptions options) {
        if (response == null || predictors == null || response.length < 2
                || predictors.length != response.length || classCount < 2
                || options == null || predictors[0] == null)
            throw new IllegalArgumentException("multinomial inputs are invalid");
        int rows = response.length, columns = predictors[0].length;
        double[] x = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            if (response[row] < 0 || response[row] >= classCount
                    || predictors[row] == null || predictors[row].length != columns)
                throw new IllegalArgumentException("response classes or predictor rows are invalid");
            for (int column = 0; column < columns; column++) {
                x[row * columns + column] = predictors[row][column];
                if (!Double.isFinite(x[row * columns + column])) throw new IllegalArgumentException("predictors must be finite");
            }
        }
        double[] beta = new double[(classCount - 1) * columns];
        double likelihood = logLikelihood(response, x, beta, rows, columns, classCount);
        boolean converged = false; int iterations = 0;
        for (int iteration = 1; iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            double[] gradient = gradient(response, x, beta, rows, columns, classCount);
            double step = options.initialStep(); double[] candidate = beta.clone();
            double candidateLikelihood = Double.NEGATIVE_INFINITY;
            for (int attempt = 0; attempt < 40; attempt++) {
                for (int i = 0; i < beta.length; i++) candidate[i] = beta[i] + step * gradient[i];
                candidateLikelihood = logLikelihood(response, x, candidate, rows, columns, classCount);
                if (candidateLikelihood >= likelihood) break;
                step *= 0.5;
            }
            if (candidateLikelihood < likelihood) {
                candidate = beta.clone();
                candidateLikelihood = likelihood;
            }
            double change = 0.0;
            for (int i = 0; i < beta.length; i++) change = Math.max(change, Math.abs(candidate[i] - beta[i]) / (1 + Math.abs(beta[i])));
            beta = candidate; double improvement = candidateLikelihood - likelihood; likelihood = candidateLikelihood;
            if (change <= options.relativeTolerance() && improvement <= options.relativeTolerance()) { converged = true; break; }
        }
        return new MultinomialRegressionResult(beta, probabilities(x, beta, rows, columns, classCount),
            likelihood, rows, columns, classCount, iterations, converged);
    }

    public static MultinomialRegressionResult fit(int[] response, double[][] predictors, int classCount) {
        return fit(response, predictors, classCount, MultinomialOptions.defaults());
    }

    private static double[] gradient(int[] y, double[] x, double[] beta, int rows, int columns, int classes) {
        double[] result = new double[beta.length], probabilities = probabilities(x, beta, rows, columns, classes);
        for (int row = 0; row < rows; row++) for (int cls = 1; cls < classes; cls++) {
            double residual = (y[row] == cls ? 1.0 : 0.0) - probabilities[row * classes + cls];
            for (int column = 0; column < columns; column++) result[(cls - 1) * columns + column] += x[row * columns + column] * residual;
        }
        return result;
    }
    private static double logLikelihood(int[] y, double[] x, double[] beta, int rows, int columns, int classes) {
        double[] probabilities = probabilities(x, beta, rows, columns, classes); double result = 0;
        for (int row = 0; row < rows; row++) result += Math.log(Math.max(1e-300, probabilities[row * classes + y[row]]));
        return result;
    }
    private static double[] probabilities(double[] x, double[] beta, int rows, int columns, int classes) {
        double[] result = new double[rows * classes];
        for (int row = 0; row < rows; row++) {
            double maximum = 0.0;
            for (int cls = 1; cls < classes; cls++) { double value = 0; for (int column = 0; column < columns; column++) value += x[row * columns + column] * beta[(cls - 1) * columns + column]; maximum = Math.max(maximum, value); }
            double denominator = Math.exp(-maximum);
            for (int cls = 1; cls < classes; cls++) { double value = 0; for (int column = 0; column < columns; column++) value += x[row * columns + column] * beta[(cls - 1) * columns + column]; result[row * classes + cls] = Math.exp(value - maximum); denominator += result[row * classes + cls]; }
            result[row * classes] = Math.exp(-maximum) / denominator;
            for (int cls = 1; cls < classes; cls++) result[row * classes + cls] /= denominator;
        }
        return result;
    }
}
