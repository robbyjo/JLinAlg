/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import jdistlib.Normal;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.MatrixOps;

/**
 * Fast maximum-likelihood beta regression with modeled mean and precision.
 *
 * <p>The implementation follows the mean/precision parameterization and
 * expected-information scoring equations used by R's {@code betareg}. Its
 * specialized two-block, row-major kernel avoids the generic distributional
 * model's per-observation parameter dispatch.</p>
 */
public final class BetaRegression {
    private BetaRegression() { }

    /** Fits a logit mean model with identity-linked constant precision. */
    public static BetaRegressionResult fit(
            double[] response, double[][] meanDesign) {
        return fit(response, meanDesign, intercept(response.length),
            BetaRegressionOptions.constantPrecisionDefaults(),
            BackendPolicy.CPU);
    }

    /** Fits mean and precision models using betareg-compatible default links. */
    public static BetaRegressionResult fit(
            double[] response,
            double[][] meanDesign,
            double[][] precisionDesign) {
        return fit(response, meanDesign, precisionDesign,
            BetaRegressionOptions.variablePrecisionDefaults(),
            BackendPolicy.CPU);
    }

    /** Fits conventional rectangular Java design matrices. */
    public static BetaRegressionResult fit(
            double[] response,
            double[][] meanDesign,
            double[][] precisionDesign,
            BetaRegressionOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || response.length == 0) {
            throw new IllegalArgumentException("response is required");
        }
        double[] mean = MatrixOps.rowMajor(meanDesign, response.length);
        double[] precision = MatrixOps.rowMajor(
            precisionDesign, response.length);
        return fit(response, mean, meanDesign[0].length,
            precision, precisionDesign[0].length, options, backendPolicy);
    }

    /** Fits contiguous row-major mean and precision design matrices. */
    public static BetaRegressionResult fit(
            double[] response,
            double[] meanDesign,
            int meanColumns,
            double[] precisionDesign,
            int precisionColumns,
            BetaRegressionOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || response.length == 0
                || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "response, controls, and backend are required");
        }
        int rows = response.length;
        MatrixOps.validateModelData(response, meanDesign, rows, meanColumns);
        if (precisionDesign == null
                || precisionDesign.length != rows * precisionColumns
                || precisionColumns < 1) {
            throw new IllegalArgumentException(
                "precision design dimensions are invalid");
        }
        MatrixOps.requireFinite(precisionDesign, "precision design");
        for (double value : response) {
            if (!(value > 0.0 && value < 1.0)) {
                throw new IllegalArgumentException(
                    "beta responses must lie strictly between zero and one");
            }
        }
        if (rows <= meanColumns + precisionColumns) {
            throw new IllegalArgumentException(
                "beta regression requires more observations than parameters");
        }

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fit(response, meanDesign, meanColumns,
                precisionDesign, precisionColumns, options,
                context.backend(), context.provenance());
        }
    }

    private static BetaRegressionResult fit(
            double[] response,
            double[] meanDesign,
            int meanColumns,
            double[] precisionDesign,
            int precisionColumns,
            BetaRegressionOptions options,
            ComputeBackend backend,
            org.jlinalg.compute.BackendProvenance provenance) {
        int rows = response.length;
        int columns = meanColumns + precisionColumns;
        double[] coefficients = initialize(response, meanDesign, meanColumns,
            precisionDesign, precisionColumns, options, backend);
        Evaluation current = evaluate(response, meanDesign, meanColumns,
            precisionDesign, precisionColumns, coefficients, options);
        boolean converged = false;
        String message = "maximum Fisher-scoring iterations reached";
        int iterations = 0;
        SystemState system = null;

        for (int iteration = 1;
                iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            system = system(response, meanDesign, meanColumns,
                precisionDesign, precisionColumns, current, options);
            double[] step;
            try {
                step = backend.dpotrf(system.information(), columns)
                    .solve(system.gradient());
            } catch (IllegalArgumentException exception) {
                message = "expected information is not positive definite";
                break;
            }
            limit(step, options.maximumStep());
            double[] candidateCoefficients = null;
            Evaluation candidate = null;
            double scale = 1.0;
            for (int attempt = 0; attempt < 40; attempt++) {
                double[] trial = coefficients.clone();
                for (int column = 0; column < columns; column++) {
                    trial[column] += scale * step[column];
                }
                Evaluation trialEvaluation = evaluate(response,
                    meanDesign, meanColumns, precisionDesign,
                    precisionColumns, trial, options);
                if (Double.isFinite(trialEvaluation.logLikelihood())
                        && trialEvaluation.logLikelihood()
                            >= current.logLikelihood() - 1e-12
                                * (1.0 + Math.abs(current.logLikelihood()))) {
                    candidateCoefficients = trial;
                    candidate = trialEvaluation;
                    break;
                }
                scale *= 0.5;
            }
            if (candidate == null) {
                message = "step halving could not improve log-likelihood";
                break;
            }
            double coefficientChange = relativeMaximumChange(
                coefficients, candidateCoefficients);
            double likelihoodChange = Math.abs(candidate.logLikelihood()
                - current.logLikelihood())
                / (1.0 + Math.abs(current.logLikelihood()));
            coefficients = candidateCoefficients;
            current = candidate;
            if (coefficientChange <= options.relativeTolerance()
                    && likelihoodChange <= options.relativeTolerance()) {
                converged = true;
                message = "coefficient and likelihood tolerances reached";
                break;
            }
        }

        system = system(response, meanDesign, meanColumns,
            precisionDesign, precisionColumns, current, options);
        CholeskyFactor factor = backend.dpotrf(system.information(), columns);
        double[] covariance = factor.solve(
            MatrixOps.identity(columns), columns);
        double[] standardErrors = new double[columns];
        double[] statistics = new double[columns];
        double[] pValues = new double[columns];
        for (int column = 0; column < columns; column++) {
            standardErrors[column] = Math.sqrt(Math.max(0.0,
                covariance[column * columns + column]));
            statistics[column] = standardErrors[column] == 0.0
                ? Math.copySign(Double.POSITIVE_INFINITY, coefficients[column])
                : coefficients[column] / standardErrors[column];
            pValues[column] = Math.min(1.0, 2.0 * Normal.cumulative(
                Math.abs(statistics[column]), 0.0, 1.0, false, false));
        }
        return new BetaRegressionResult(
            Arrays.copyOfRange(coefficients, 0, meanColumns),
            Arrays.copyOfRange(coefficients, meanColumns, columns),
            covariance, standardErrors, statistics, pValues,
            current.means(), current.precisions(), current.logLikelihood(),
            rows, iterations, converged, message,
            options.meanLink(), options.precisionLink(), provenance);
    }

    private static double[] initialize(
            double[] response,
            double[] meanDesign,
            int meanColumns,
            double[] precisionDesign,
            int precisionColumns,
            BetaRegressionOptions options,
            ComputeBackend backend) {
        int rows = response.length;
        double[] target = new double[rows];
        for (int row = 0; row < rows; row++) {
            target[row] = options.meanLink().link(response[row]);
        }
        double[] meanCoefficients = LeastSquaresSolver.solve(
            meanDesign, target, rows, meanColumns, false, backend)
            .coefficients();
        double[] coefficients = new double[meanColumns + precisionColumns];
        System.arraycopy(meanCoefficients, 0, coefficients, 0, meanColumns);

        double mean = Arrays.stream(response).average().orElseThrow();
        double variance = 0.0;
        for (double value : response) {
            double difference = value - mean;
            variance += difference * difference;
        }
        variance /= Math.max(1, rows - 1);
        double precision = variance == 0.0 ? 100.0
            : Math.max(1.0, mean * (1.0 - mean) / variance - 1.0);
        double first = precisionDesign[0];
        boolean constant = first != 0.0;
        for (int row = 1; row < rows; row++) {
            constant &= precisionDesign[row * precisionColumns] == first;
        }
        if (constant) {
            coefficients[meanColumns] = options.precisionLink().link(precision)
                / first;
        }
        return coefficients;
    }

    private static Evaluation evaluate(
            double[] response,
            double[] meanDesign,
            int meanColumns,
            double[] precisionDesign,
            int precisionColumns,
            double[] coefficients,
            BetaRegressionOptions options) {
        int rows = response.length;
        double[] means = new double[rows];
        double[] precisions = new double[rows];
        double[] meanPredictors = new double[rows];
        double[] precisionPredictors = new double[rows];
        double logLikelihood = 0.0;
        for (int row = 0; row < rows; row++) {
            int meanOffset = row * meanColumns;
            double meanPredictor = 0.0;
            for (int column = 0; column < meanColumns; column++) {
                meanPredictor += meanDesign[meanOffset + column]
                    * coefficients[column];
            }
            int precisionOffset = row * precisionColumns;
            double precisionPredictor = 0.0;
            for (int column = 0; column < precisionColumns; column++) {
                precisionPredictor += precisionDesign[precisionOffset + column]
                    * coefficients[meanColumns + column];
            }
            double mean = options.meanLink().inverse(meanPredictor);
            double precision = options.precisionLink().inverse(
                precisionPredictor);
            double alpha = mean * precision;
            double beta = (1.0 - mean) * precision;
            double contribution = SpecialFunctions.logGamma(precision)
                - SpecialFunctions.logGamma(alpha)
                - SpecialFunctions.logGamma(beta)
                + (alpha - 1.0) * Math.log(response[row])
                + (beta - 1.0) * Math.log1p(-response[row]);
            if (!Double.isFinite(contribution)) {
                return new Evaluation(meanPredictors, precisionPredictors,
                    means, precisions, Double.NaN);
            }
            meanPredictors[row] = meanPredictor;
            precisionPredictors[row] = precisionPredictor;
            means[row] = mean;
            precisions[row] = precision;
            logLikelihood += contribution;
        }
        return new Evaluation(meanPredictors, precisionPredictors,
            means, precisions, logLikelihood);
    }

    private static SystemState system(
            double[] response,
            double[] meanDesign,
            int meanColumns,
            double[] precisionDesign,
            int precisionColumns,
            Evaluation evaluation,
            BetaRegressionOptions options) {
        int columns = meanColumns + precisionColumns;
        double[] gradient = new double[columns];
        double[] information = new double[columns * columns];
        for (int row = 0; row < response.length; row++) {
            double mean = evaluation.means()[row];
            double precision = evaluation.precisions()[row];
            double alpha = mean * precision;
            double beta = (1.0 - mean) * precision;
            double digammaAlpha = SpecialFunctions.digamma(alpha);
            double digammaBeta = SpecialFunctions.digamma(beta);
            double transformedResponse = Math.log(response[row])
                - Math.log1p(-response[row]);
            double transformedMean = digammaAlpha - digammaBeta;
            double meanDerivative = options.meanLink().derivative(
                evaluation.meanPredictors()[row], mean);
            double precisionDerivative = options.precisionLink().derivative(
                evaluation.precisionPredictors()[row], precision);
            double meanScore = precision
                * (transformedResponse - transformedMean) * meanDerivative;
            double precisionScore = (mean
                * (transformedResponse - transformedMean)
                + Math.log1p(-response[row]) - digammaBeta
                + SpecialFunctions.digamma(precision)) * precisionDerivative;

            double trigammaAlpha = SpecialFunctions.trigamma(alpha);
            double trigammaBeta = SpecialFunctions.trigamma(beta);
            double sumTrigamma = trigammaAlpha + trigammaBeta;
            double meanWeight = precision * precision * sumTrigamma
                * meanDerivative * meanDerivative;
            double precisionWeight = (mean * mean * trigammaAlpha
                + (1.0 - mean) * (1.0 - mean) * trigammaBeta
                - SpecialFunctions.trigamma(precision))
                * precisionDerivative * precisionDerivative;
            double crossWeight = precision
                * (mean * sumTrigamma - trigammaBeta)
                * meanDerivative * precisionDerivative;

            int meanOffset = row * meanColumns;
            for (int left = 0; left < meanColumns; left++) {
                double leftValue = meanDesign[meanOffset + left];
                gradient[left] += leftValue * meanScore;
                for (int right = 0; right <= left; right++) {
                    information[left * columns + right] += leftValue
                        * meanWeight * meanDesign[meanOffset + right];
                }
                int precisionOffset = row * precisionColumns;
                for (int right = 0; right < precisionColumns; right++) {
                    information[left * columns + meanColumns + right] +=
                        leftValue * crossWeight
                            * precisionDesign[precisionOffset + right];
                }
            }
            int precisionOffset = row * precisionColumns;
            for (int left = 0; left < precisionColumns; left++) {
                int leftIndex = meanColumns + left;
                double leftValue = precisionDesign[precisionOffset + left];
                gradient[leftIndex] += leftValue * precisionScore;
                for (int right = 0; right <= left; right++) {
                    int rightIndex = meanColumns + right;
                    information[leftIndex * columns + rightIndex] += leftValue
                        * precisionWeight
                        * precisionDesign[precisionOffset + right];
                }
            }
        }
        for (int row = 0; row < columns; row++) {
            for (int column = 0; column < row; column++) {
                if (row >= meanColumns && column < meanColumns) {
                    information[row * columns + column] =
                        information[column * columns + row];
                } else {
                    information[column * columns + row] =
                        information[row * columns + column];
                }
            }
        }
        return new SystemState(gradient, information);
    }

    private static double relativeMaximumChange(
            double[] current, double[] candidate) {
        double maximum = 0.0;
        for (int index = 0; index < current.length; index++) {
            maximum = Math.max(maximum,
                Math.abs(candidate[index] - current[index])
                    / (1.0 + Math.abs(current[index])));
        }
        return maximum;
    }

    private static void limit(double[] values, double maximum) {
        double observed = 0.0;
        for (double value : values) {
            observed = Math.max(observed, Math.abs(value));
        }
        if (observed > maximum) {
            double scale = maximum / observed;
            for (int index = 0; index < values.length; index++) {
                values[index] *= scale;
            }
        }
    }

    private static double[][] intercept(int rows) {
        double[][] result = new double[rows][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }

    private record Evaluation(
            double[] meanPredictors,
            double[] precisionPredictors,
            double[] means,
            double[] precisions,
            double logLikelihood) { }

    private record SystemState(double[] gradient, double[] information) { }
}
