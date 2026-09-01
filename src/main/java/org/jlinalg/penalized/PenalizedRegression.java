/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlinalg.internal.MatrixOps;

/** Fast Gaussian ridge, LASSO, and elastic-net regression. */
public final class PenalizedRegression {
    private static final double RIDGE_PATH_ALPHA = 1e-3;

    private PenalizedRegression() { }

    /** Fits ridge regression ({@code alpha = 0}). */
    public static PenalizedRegressionResult ridge(
            double[] response, double[][] predictors, double lambda) {
        return ridge(response, predictors, lambda, ElasticNetOptions.defaults());
    }

    /** Fits ridge regression while preserving all non-alpha options. */
    public static PenalizedRegressionResult ridge(
            double[] response,
            double[][] predictors,
            double lambda,
            ElasticNetOptions options) {
        return fit(response, predictors, lambda, copyWithAlpha(options, 0.0));
    }

    /** Fits LASSO regression ({@code alpha = 1}). */
    public static PenalizedRegressionResult lasso(
            double[] response, double[][] predictors, double lambda) {
        return lasso(response, predictors, lambda, ElasticNetOptions.defaults());
    }

    /** Fits LASSO while preserving all non-alpha options. */
    public static PenalizedRegressionResult lasso(
            double[] response,
            double[][] predictors,
            double lambda,
            ElasticNetOptions options) {
        return fit(response, predictors, lambda, copyWithAlpha(options, 1.0));
    }

    /** Fits an elastic net using {@link ElasticNetOptions#alpha()}. */
    public static PenalizedRegressionResult fit(
            double[] response,
            double[][] predictors,
            double lambda,
            ElasticNetOptions options) {
        return path(response, predictors, new double[] {lambda}, options).fit(0);
    }

    /** Fits a descending lambda sequence with warm starts. */
    public static PenalizedRegressionPath path(
            double[] response,
            double[][] predictors,
            double[] lambdas,
            ElasticNetOptions options) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] design = MatrixOps.rowMajor(predictors, response.length);
        return path(response, design, response.length, predictors[0].length,
            lambdas, options);
    }

    /** Fits a path from a contiguous row-major predictor matrix. */
    public static PenalizedRegressionPath path(
            double[] response,
            double[] predictors,
            int rows,
            int columns,
            double[] lambdas,
            ElasticNetOptions options) {
        MatrixOps.validateModelData(response, predictors, rows, columns);
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        double[] penalties = validateLambdas(lambdas);
        PreparedData data = prepare(
            response, predictors, rows, columns, options);
        return fitPrepared(data, penalties, options);
    }

    /**
     * Generates a log-linear lambda path from the all-zero threshold down to
     * {@code minimumRatio * lambdaMax} and fits it with warm starts.
     */
    public static PenalizedRegressionPath automaticPath(
            double[] response,
            double[][] predictors,
            int lambdaCount,
            double minimumRatio,
            ElasticNetOptions options) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        if (lambdaCount < 2) {
            throw new IllegalArgumentException("lambdaCount must be at least two");
        }
        if (!Double.isFinite(minimumRatio)
                || !(minimumRatio > 0.0) || !(minimumRatio < 1.0)) {
            throw new IllegalArgumentException("minimumRatio must lie in (0, 1)");
        }
        double[] design = MatrixOps.rowMajor(predictors, response.length);
        int columns = predictors[0].length;
        MatrixOps.validateModelData(
            response, design, response.length, columns);
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        PreparedData data = prepare(
            response, design, response.length, columns, options);
        double maximum = lambdaMaximum(data, options.alpha());
        if (!(maximum > 0.0) || !Double.isFinite(maximum)) {
            maximum = 1.0;
        }
        double[] lambdas = new double[lambdaCount];
        double logMaximum = Math.log(maximum);
        double logMinimum = Math.log(maximum * minimumRatio);
        for (int index = 0; index < lambdaCount; index++) {
            double fraction = (double) index / (lambdaCount - 1);
            lambdas[index] = Math.exp(
                logMaximum + fraction * (logMinimum - logMaximum));
        }
        return fitPrepared(data, lambdas, options);
    }

    static PenalizedRegressionPath fitPrepared(
            PreparedData data,
            double[] lambdas,
            ElasticNetOptions options) {
        int rows = data.rows();
        int columns = data.columns();
        double[] beta = new double[columns];
        double[] residual = data.workingResponse().clone();
        List<PenalizedRegressionResult> results = new ArrayList<>(lambdas.length);

        for (double lambda : lambdas) {
            boolean converged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumIterations(); iteration++) {
                iterations = iteration;
                double maximumChange = 0.0;
                double maximumCoefficient = 0.0;
                for (int column = 0; column < columns; column++) {
                    double old = beta[column];
                    if (old != 0.0) {
                        for (int row = 0; row < rows; row++) {
                            residual[row] += data.workingPredictors()[
                                row * columns + column] * old;
                        }
                    }
                    double correlation = 0.0;
                    for (int row = 0; row < rows; row++) {
                        correlation += data.weights()[row]
                            * data.workingPredictors()[row * columns + column]
                            * residual[row];
                    }
                    correlation /= rows;
                    double penaltyFactor = data.penaltyFactors()[column];
                    double numerator = softThreshold(
                        correlation, lambda * options.alpha() * penaltyFactor);
                    double denominator = data.columnSquaredNorms()[column]
                        + lambda * (1.0 - options.alpha()) * penaltyFactor;
                    double updated = denominator > 0.0
                        ? numerator / denominator : 0.0;
                    beta[column] = updated;
                    if (updated != 0.0) {
                        for (int row = 0; row < rows; row++) {
                            residual[row] -= data.workingPredictors()[
                                row * columns + column] * updated;
                        }
                    }
                    maximumChange = Math.max(
                        maximumChange, Math.abs(updated - old));
                    maximumCoefficient = Math.max(
                        maximumCoefficient, Math.abs(updated));
                }
                if (maximumChange <= options.relativeTolerance()
                        * (1.0 + maximumCoefficient)) {
                    converged = true;
                    break;
                }
            }
            results.add(result(data, beta, lambda, options.alpha(),
                iterations, converged));
        }
        return new PenalizedRegressionPath(results);
    }

    static PreparedData prepare(
            double[] response,
            double[] predictors,
            int rows,
            int columns,
            ElasticNetOptions options) {
        double[] weights = prepareWeights(options.observationWeights(), rows);
        double[] penaltyFactors = preparePenaltyFactors(
            options.penaltyFactors(), columns);
        double responseMean = options.fitIntercept()
            ? weightedMean(response, weights, rows) : 0.0;
        double[] predictorMeans = new double[columns];
        if (options.fitIntercept()) {
            for (int column = 0; column < columns; column++) {
                double sum = 0.0;
                for (int row = 0; row < rows; row++) {
                    sum += weights[row] * predictors[row * columns + column];
                }
                predictorMeans[column] = sum / rows;
            }
        }

        double[] scales = new double[columns];
        Arrays.fill(scales, 1.0);
        double[] workingPredictors = new double[predictors.length];
        for (int column = 0; column < columns; column++) {
            if (options.standardize()) {
                double sumSquares = 0.0;
                for (int row = 0; row < rows; row++) {
                    double centered = predictors[row * columns + column]
                        - predictorMeans[column];
                    sumSquares += weights[row] * centered * centered;
                }
                double scale = Math.sqrt(sumSquares / rows);
                scales[column] = scale > 1e-14 ? scale : 1.0;
            }
            for (int row = 0; row < rows; row++) {
                workingPredictors[row * columns + column] =
                    (predictors[row * columns + column] - predictorMeans[column])
                        / scales[column];
            }
        }
        double[] workingResponse = new double[rows];
        for (int row = 0; row < rows; row++) {
            workingResponse[row] = response[row] - responseMean;
        }
        double[] squaredNorms = new double[columns];
        for (int column = 0; column < columns; column++) {
            for (int row = 0; row < rows; row++) {
                double value = workingPredictors[row * columns + column];
                squaredNorms[column] += weights[row] * value * value / rows;
            }
        }
        return new PreparedData(response.clone(), predictors.clone(),
            workingResponse, workingPredictors, weights, penaltyFactors,
            predictorMeans, scales, squaredNorms, responseMean, rows, columns);
    }

    static double lambdaMaximum(PreparedData data, double alpha) {
        double effectiveAlpha = Math.max(alpha, RIDGE_PATH_ALPHA);
        double maximum = 0.0;
        for (int column = 0; column < data.columns(); column++) {
            double correlation = 0.0;
            for (int row = 0; row < data.rows(); row++) {
                correlation += data.weights()[row]
                    * data.workingPredictors()[row * data.columns() + column]
                    * data.workingResponse()[row];
            }
            correlation = Math.abs(correlation / data.rows());
            maximum = Math.max(maximum,
                correlation / (effectiveAlpha * data.penaltyFactors()[column]));
        }
        return maximum;
    }

    private static PenalizedRegressionResult result(
            PreparedData data,
            double[] workingBeta,
            double lambda,
            double alpha,
            int iterations,
            boolean converged) {
        double[] coefficients = new double[data.columns()];
        double intercept = data.responseMean();
        for (int column = 0; column < data.columns(); column++) {
            coefficients[column] = workingBeta[column] / data.scales()[column];
            intercept -= data.predictorMeans()[column] * coefficients[column];
        }
        double[] fitted = new double[data.rows()];
        double[] residuals = new double[data.rows()];
        double weightedRss = 0.0;
        for (int row = 0; row < data.rows(); row++) {
            double value = intercept;
            for (int column = 0; column < data.columns(); column++) {
                value += data.originalPredictors()[row * data.columns() + column]
                    * coefficients[column];
            }
            fitted[row] = value;
            residuals[row] = data.originalResponse()[row] - value;
            weightedRss += data.weights()[row]
                * residuals[row] * residuals[row];
        }
        double l1 = 0.0;
        double l2 = 0.0;
        int active = 0;
        for (int column = 0; column < data.columns(); column++) {
            double value = workingBeta[column];
            l1 += data.penaltyFactors()[column] * Math.abs(value);
            l2 += data.penaltyFactors()[column] * value * value;
            if (value != 0.0) {
                active++;
            }
        }
        double objective = weightedRss / (2.0 * data.rows())
            + lambda * (alpha * l1 + 0.5 * (1.0 - alpha) * l2);
        return new PenalizedRegressionResult(
            intercept, coefficients, fitted, residuals,
            lambda, alpha, objective, weightedRss, active,
            iterations, converged,
            converged ? "coordinate tolerance reached"
                : "maximum coordinate-descent iterations reached");
    }

    private static double[] prepareWeights(double[] supplied, int rows) {
        double[] weights;
        if (supplied == null) {
            weights = new double[rows];
            Arrays.fill(weights, 1.0);
        } else {
            if (supplied.length != rows) {
                throw new IllegalArgumentException(
                    "observation weight length must equal rows");
            }
            weights = supplied.clone();
        }
        double sum = 0.0;
        for (double value : weights) {
            if (!Double.isFinite(value) || !(value > 0.0)) {
                throw new IllegalArgumentException(
                    "observation weights must be finite and positive");
            }
            sum += value;
        }
        double scale = rows / sum;
        for (int index = 0; index < rows; index++) {
            weights[index] *= scale;
        }
        return weights;
    }

    private static double[] preparePenaltyFactors(
            double[] supplied, int columns) {
        double[] factors;
        if (supplied == null) {
            factors = new double[columns];
            Arrays.fill(factors, 1.0);
        } else {
            if (supplied.length != columns) {
                throw new IllegalArgumentException(
                    "penalty factor length must equal predictor columns");
            }
            factors = supplied.clone();
        }
        for (double value : factors) {
            if (!Double.isFinite(value) || !(value > 0.0)) {
                throw new IllegalArgumentException(
                    "penalty factors must be finite and positive");
            }
        }
        return factors;
    }

    private static double weightedMean(
            double[] values, double[] weights, int rows) {
        double sum = 0.0;
        for (int row = 0; row < rows; row++) {
            sum += weights[row] * values[row];
        }
        return sum / rows;
    }

    private static double softThreshold(double value, double threshold) {
        if (value > threshold) {
            return value - threshold;
        }
        if (value < -threshold) {
            return value + threshold;
        }
        return 0.0;
    }

    private static double[] validateLambdas(double[] lambdas) {
        if (lambdas == null || lambdas.length == 0) {
            throw new IllegalArgumentException("at least one lambda is required");
        }
        double[] result = lambdas.clone();
        for (int index = 0; index < result.length; index++) {
            if (!Double.isFinite(result[index]) || result[index] < 0.0) {
                throw new IllegalArgumentException(
                    "lambdas must be finite and nonnegative");
            }
            if (index > 0 && !(result[index] < result[index - 1])) {
                throw new IllegalArgumentException(
                    "lambdas must be in strictly descending order");
            }
        }
        return result;
    }

    static ElasticNetOptions copyWithAlpha(
            ElasticNetOptions source, double alpha) {
        if (source == null) {
            throw new IllegalArgumentException("options are required");
        }
        ElasticNetOptions.Builder builder = ElasticNetOptions.builder()
            .alpha(alpha)
            .fitIntercept(source.fitIntercept())
            .standardize(source.standardize())
            .maximumIterations(source.maximumIterations())
            .relativeTolerance(source.relativeTolerance());
        if (source.observationWeights() != null) {
            builder.observationWeights(source.observationWeights());
        }
        if (source.penaltyFactors() != null) {
            builder.penaltyFactors(source.penaltyFactors());
        }
        return builder.build();
    }

    record PreparedData(
            double[] originalResponse,
            double[] originalPredictors,
            double[] workingResponse,
            double[] workingPredictors,
            double[] weights,
            double[] penaltyFactors,
            double[] predictorMeans,
            double[] scales,
            double[] columnSquaredNorms,
            double responseMean,
            int rows,
            int columns) {
    }
}
