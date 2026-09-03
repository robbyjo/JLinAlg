/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import org.jlinalg.internal.MatrixOps;

/** Fast Gaussian ridge, LASSO, and elastic-net regression. */
public final class PenalizedRegression {
    private static final double RIDGE_PATH_ALPHA = 1e-3;
    private static final int COVARIANCE_MAX_COLUMNS = 2_048;

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
        PreparedData data = prepareData(
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
        double[] design = MatrixOps.rowMajor(predictors, response.length);
        return automaticPath(response, design, response.length,
            predictors[0].length, lambdaCount, minimumRatio, options);
    }

    /** Generates and fits a path from a contiguous row-major matrix. */
    public static PenalizedRegressionPath automaticPath(
            double[] response,
            double[] predictors,
            int rows,
            int columns,
            int lambdaCount,
            double minimumRatio,
            ElasticNetOptions options) {
        MatrixOps.validateModelData(response, predictors, rows, columns);
        if (lambdaCount < 2) {
            throw new IllegalArgumentException("lambdaCount must be at least two");
        }
        if (!Double.isFinite(minimumRatio)
                || !(minimumRatio > 0.0) || !(minimumRatio < 1.0)) {
            throw new IllegalArgumentException("minimumRatio must lie in (0, 1)");
        }
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        PreparedData data = prepareData(
            response, predictors, rows, columns, options);
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
        double[] residual = data.gramMatrix() == null
            ? data.workingResponse().clone() : null;
        double[] correlations = data.gramMatrix() == null
            ? null : data.responseCorrelations().clone();
        List<PenalizedRegressionResult> results = new ArrayList<>(lambdas.length);

        for (double lambda : lambdas) {
            boolean converged = false;
            int iterations = 0;
            for (int iteration = 1;
                    iteration <= options.maximumIterations(); iteration++) {
                iterations = iteration;
                double maximumContribution = 0.0;
                for (int column = 0; column < columns; column++) {
                    double old = beta[column];
                    double correlation;
                    if (correlations != null) {
                        correlation = correlations[column]
                            + data.columnSquaredNorms()[column] * old;
                    } else {
                        int offset = column * rows;
                        if (old != 0.0) {
                            for (int row = 0; row < rows; row++) {
                                residual[row] += data.workingPredictors()[
                                    offset + row] * old;
                            }
                        }
                        correlation = 0.0;
                        for (int row = 0; row < rows; row++) {
                            correlation += data.weights()[row]
                                * data.workingPredictors()[offset + row]
                                * residual[row];
                        }
                        correlation /= rows;
                    }
                    double penaltyFactor = data.penaltyFactors()[column];
                    double numerator = softThreshold(
                        correlation, lambda * options.alpha() * penaltyFactor);
                    double denominator = data.columnSquaredNorms()[column]
                        + lambda * (1.0 - options.alpha()) * penaltyFactor;
                    double updated = denominator > 0.0
                        ? numerator / denominator : 0.0;
                    beta[column] = updated;
                    if (correlations != null) {
                        double change = updated - old;
                        if (change != 0.0) {
                            for (int target = 0; target < columns; target++) {
                                correlations[target] -= data.gramMatrix()[
                                    target * columns + column] * change;
                            }
                        }
                    } else if (updated != 0.0) {
                        int offset = column * rows;
                        for (int row = 0; row < rows; row++) {
                            residual[row] -= data.workingPredictors()[
                                offset + row] * updated;
                        }
                    }
                    double change = updated - old;
                    maximumContribution = Math.max(maximumContribution,
                        data.columnSquaredNorms()[column] * change * change);
                }
                if (maximumContribution <= options.relativeTolerance()
                        * data.nullDeviance()) {
                    converged = true;
                    break;
                }
            }
            results.add(result(data, beta, lambda, options.alpha(),
                iterations, converged));
        }
        return new PenalizedRegressionPath(results);
    }

    /** Preprocesses a conventional matrix for reuse across penalty paths. */
    public static Prepared prepare(
            double[] response,
            double[][] predictors,
            ElasticNetOptions options) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] design = MatrixOps.rowMajor(predictors, response.length);
        return prepare(response, design, response.length,
            predictors[0].length, options);
    }

    /** Preprocesses a contiguous row-major matrix for reusable path fits. */
    public static Prepared prepare(
            double[] response,
            double[] predictors,
            int rows,
            int columns,
            ElasticNetOptions options) {
        MatrixOps.validateModelData(response, predictors, rows, columns);
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        return new Prepared(prepareData(
            response, predictors, rows, columns, options), options);
    }

    static PreparedData prepareData(
            double[] response,
            double[] predictors,
            int rows,
            int columns,
            ElasticNetOptions options) {
        double[] weights = prepareWeights(options.observationWeights(), rows);
        boolean unitWeights = unitWeights(weights);
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
        double[] squaredNorms = new double[columns];
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
            int offset = column * rows;
            for (int row = 0; row < rows; row++) {
                double value =
                    (predictors[row * columns + column] - predictorMeans[column])
                        / scales[column];
                workingPredictors[offset + row] = value;
                squaredNorms[column] += (unitWeights ? 1.0 : weights[row])
                    * value * value / rows;
            }
        }
        double[] workingResponse = new double[rows];
        double nullDeviance = 0.0;
        for (int row = 0; row < rows; row++) {
            workingResponse[row] = response[row] - responseMean;
            nullDeviance += weights[row]
                * workingResponse[row] * workingResponse[row] / rows;
        }
        double[] responseCorrelations = new double[columns];
        for (int column = 0; column < columns; column++) {
            int offset = column * rows;
            for (int row = 0; row < rows; row++) {
                responseCorrelations[column] += (unitWeights ? 1.0 : weights[row])
                    * workingPredictors[offset + row] * workingResponse[row]
                    / rows;
            }
        }
        double[] gram = columns <= rows && columns <= COVARIANCE_MAX_COLUMNS
            ? gram(workingPredictors, weights, unitWeights, rows, columns,
                options.parallelism()) : null;
        return new PreparedData(response.clone(), predictors.clone(),
            workingResponse, workingPredictors, weights, penaltyFactors,
            predictorMeans, scales, squaredNorms, responseCorrelations, gram,
            responseMean, nullDeviance, rows, columns);
    }

    private static double[] gram(double[] predictors, double[] weights,
            boolean unitWeights, int rows, int columns, int parallelism) {
        double[] result = new double[columns * columns];
        Runnable compute = () -> IntStream.range(0, columns).parallel()
            .forEach(column -> {
                int columnOffset = column * rows;
                for (int other = 0; other <= column; other++) {
                    int otherOffset = other * rows;
                    double value = 0.0;
                    if (unitWeights) {
                        for (int row = 0; row < rows; row++) {
                            value += predictors[columnOffset + row]
                                * predictors[otherOffset + row];
                        }
                    } else {
                        for (int row = 0; row < rows; row++) {
                            value += weights[row]
                                * predictors[columnOffset + row]
                                * predictors[otherOffset + row];
                        }
                    }
                    value /= rows;
                    result[column * columns + other] = value;
                    result[other * columns + column] = value;
                }
            });
        if (parallelism == 1) {
            for (int column = 0; column < columns; column++) {
                int columnOffset = column * rows;
                for (int other = 0; other <= column; other++) {
                    int otherOffset = other * rows;
                    double value = 0.0;
                    if (unitWeights) {
                        for (int row = 0; row < rows; row++) {
                            value += predictors[columnOffset + row]
                                * predictors[otherOffset + row];
                        }
                    } else {
                        for (int row = 0; row < rows; row++) {
                            value += weights[row]
                                * predictors[columnOffset + row]
                                * predictors[otherOffset + row];
                        }
                    }
                    value /= rows;
                    result[column * columns + other] = value;
                    result[other * columns + column] = value;
                }
            }
            return result;
        }
        ForkJoinPool pool = new ForkJoinPool(parallelism);
        try {
            pool.submit(compute).join();
        } finally {
            pool.shutdown();
        }
        return result;
    }

    static double lambdaMaximum(PreparedData data, double alpha) {
        double effectiveAlpha = Math.max(alpha, RIDGE_PATH_ALPHA);
        double[] correlations = data.responseCorrelations().clone();
        if (data.gramMatrix() != null) {
            fitUnpenalized(data, correlations);
        } else {
            correlations = residualizedCorrelations(data);
        }
        double maximum = 0.0;
        for (int column = 0; column < data.columns(); column++) {
            if (data.penaltyFactors()[column] == 0.0) {
                continue;
            }
            double correlation = Math.abs(correlations[column]);
            maximum = Math.max(maximum,
                correlation / (effectiveAlpha * data.penaltyFactors()[column]));
        }
        return maximum;
    }

    private static double[] residualizedCorrelations(PreparedData data) {
        int rows = data.rows();
        int columns = data.columns();
        double[] residual = data.workingResponse().clone();
        double[] beta = new double[columns];
        for (int iteration = 0; iteration < 100_000; iteration++) {
            double maximumChange = 0.0;
            for (int column = 0; column < columns; column++) {
                if (data.penaltyFactors()[column] != 0.0) continue;
                int offset = column * rows;
                double old = beta[column];
                if (old != 0.0) {
                    for (int row = 0; row < rows; row++) {
                        residual[row] += data.workingPredictors()[offset + row]
                            * old;
                    }
                }
                double correlation = 0.0;
                for (int row = 0; row < rows; row++) {
                    correlation += data.weights()[row]
                        * data.workingPredictors()[offset + row] * residual[row];
                }
                correlation /= rows;
                double norm = data.columnSquaredNorms()[column];
                double updated = norm > 0.0 ? correlation / norm : 0.0;
                beta[column] = updated;
                if (updated != 0.0) {
                    for (int row = 0; row < rows; row++) {
                        residual[row] -= data.workingPredictors()[offset + row]
                            * updated;
                    }
                }
                maximumChange = Math.max(
                    maximumChange, Math.abs(updated - old));
            }
            if (maximumChange <= 1e-12) break;
        }
        double[] result = new double[columns];
        for (int column = 0; column < columns; column++) {
            int offset = column * rows;
            for (int row = 0; row < rows; row++) {
                result[column] += data.weights()[row]
                    * data.workingPredictors()[offset + row] * residual[row]
                    / rows;
            }
        }
        return result;
    }

    private static void fitUnpenalized(
            PreparedData data, double[] correlations) {
        int columns = data.columns();
        double[] beta = new double[columns];
        for (int iteration = 0; iteration < 100_000; iteration++) {
            double maximumChange = 0.0;
            double maximumCoefficient = 0.0;
            for (int column = 0; column < columns; column++) {
                if (data.penaltyFactors()[column] != 0.0) continue;
                double norm = data.columnSquaredNorms()[column];
                double old = beta[column];
                double updated = norm > 0.0
                    ? (correlations[column] + norm * old) / norm : 0.0;
                double change = updated - old;
                beta[column] = updated;
                if (change != 0.0) {
                    for (int target = 0; target < columns; target++) {
                        correlations[target] -= data.gramMatrix()[
                            target * columns + column] * change;
                    }
                }
                maximumChange = Math.max(maximumChange, Math.abs(change));
                maximumCoefficient = Math.max(
                    maximumCoefficient, Math.abs(updated));
            }
            if (maximumChange <= 1e-12 * (1.0 + maximumCoefficient)) return;
        }
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
        double weightedRss = weightedRss(data, workingBeta,
            intercept, coefficients);
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
            intercept, coefficients, data.originalResponse(),
            data.originalPredictors(), data.rows(), data.columns(),
            lambda, alpha, objective, weightedRss, active,
            iterations, converged,
            converged ? "coordinate tolerance reached"
                : "maximum coordinate-descent iterations reached");
    }

    private static double weightedRss(PreparedData data, double[] beta,
            double intercept, double[] coefficients) {
        if (data.gramMatrix() != null) {
            double value = data.nullDeviance();
            for (int row = 0; row < data.columns(); row++) {
                value -= 2.0 * beta[row] * data.responseCorrelations()[row];
                for (int column = 0; column < data.columns(); column++) {
                    value += beta[row] * data.gramMatrix()[
                        row * data.columns() + column] * beta[column];
                }
            }
            return Math.max(0.0, value * data.rows());
        }
        double result = 0.0;
        for (int row = 0; row < data.rows(); row++) {
            double fitted = intercept;
            for (int column = 0; column < data.columns(); column++) {
                fitted += data.originalPredictors()[
                    row * data.columns() + column] * coefficients[column];
            }
            double residual = data.originalResponse()[row] - fitted;
            result += data.weights()[row] * residual * residual;
        }
        return result;
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
            if (!Double.isFinite(value) || value < 0.0) {
                throw new IllegalArgumentException(
                    "penalty factors must be finite and nonnegative");
            }
        }
        return factors;
    }

    private static boolean unitWeights(double[] weights) {
        for (double weight : weights) {
            if (weight != 1.0) return false;
        }
        return true;
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
            .parallelism(source.parallelism())
            .relativeTolerance(source.relativeTolerance());
        if (source.observationWeights() != null) {
            builder.observationWeights(source.observationWeights());
        }
        if (source.penaltyFactors() != null) {
            builder.penaltyFactors(source.penaltyFactors());
        }
        return builder.build();
    }

    /** Standardized design and covariance data reusable across alpha values. */
    public static final class Prepared {
        private final PreparedData data;
        private final ElasticNetOptions options;

        private Prepared(PreparedData data, ElasticNetOptions options) {
            this.data = data;
            this.options = options;
        }

        public int rows() { return data.rows(); }
        public int columns() { return data.columns(); }

        /** Fits descending penalties at the requested elastic-net alpha. */
        public PenalizedRegressionPath path(double[] lambdas, double alpha) {
            return fitPrepared(data, validateLambdas(lambdas),
                copyWithAlpha(options, alpha));
        }

        /** Generates and fits a log-linear path at the requested alpha. */
        public PenalizedRegressionPath automaticPath(
                int lambdaCount, double minimumRatio, double alpha) {
            if (lambdaCount < 2) {
                throw new IllegalArgumentException(
                    "lambdaCount must be at least two");
            }
            if (!Double.isFinite(minimumRatio)
                    || !(minimumRatio > 0.0) || !(minimumRatio < 1.0)) {
                throw new IllegalArgumentException(
                    "minimumRatio must lie in (0, 1)");
            }
            ElasticNetOptions fitOptions = copyWithAlpha(options, alpha);
            double maximum = lambdaMaximum(data, alpha);
            if (!(maximum > 0.0) || !Double.isFinite(maximum)) maximum = 1.0;
            double[] lambdas = new double[lambdaCount];
            double logMaximum = Math.log(maximum);
            double logMinimum = Math.log(maximum * minimumRatio);
            for (int index = 0; index < lambdaCount; index++) {
                double fraction = (double) index / (lambdaCount - 1);
                lambdas[index] = Math.exp(logMaximum
                    + fraction * (logMinimum - logMaximum));
            }
            return fitPrepared(data, lambdas, fitOptions);
        }
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
            double[] responseCorrelations,
            double[] gramMatrix,
            double responseMean,
            double nullDeviance,
            int rows,
            int columns) {
    }
}
