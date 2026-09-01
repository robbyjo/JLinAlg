/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

import java.util.Arrays;
import jdistlib.Normal;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.LeastSquaresSolver.Solution;
import org.jlinalg.internal.CompleteCases;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.ols.RankDeficiencyStrategy;

/** Generalized linear models fitted by stabilized iteratively reweighted least squares. */
public final class Glm {
    private static final double MINIMUM_WORKING_WEIGHT = 1e-15;
    private static final double MAXIMUM_WORKING_WEIGHT = 1e150;

    private Glm() { }

    /** Fits a GLM with unit prior weights, zero offset, and default controls. */
    public static GlmResult fit(
            double[] response, double[][] design, GlmFamily family) {
        return fit(response, design, family, null, null,
            GlmOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a GLM from a conventional rectangular Java matrix. */
    public static GlmResult fit(
            double[] response,
            double[][] design,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajorUnchecked(design, response.length);
        return fit(response, rowMajor, response.length, design[0].length,
            family, priorWeights, offset, options, backendPolicy);
    }

    /** Fits a GLM from a contiguous row-major design matrix. */
    public static GlmResult fit(
            double[] response,
            double[] design,
            int rows,
            int columns,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmOptions options,
            BackendPolicy backendPolicy) {
        if (family == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "family, options, and backendPolicy are required");
        }
        CompleteCases.Data complete = CompleteCases.prepare(
            response, design, rows, columns, priorWeights, offset,
            options.missingDataPolicy());
        int effectiveRows = complete.response().length;
        double[] weights = prepareWeights(complete.weights(), effectiveRows);
        double[] offsets = prepareOffset(complete.offset(), effectiveRows);
        validateResponses(complete.response(), weights, family);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fit(complete.response(), complete.design(),
                effectiveRows, columns, family,
                weights, offsets, options,
                complete.retainedRows(), complete.originalRows(),
                context.backend(), context.provenance());
        }
    }

    private static GlmResult fit(
            double[] response,
            double[] design,
            int rows,
            int columns,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            GlmOptions options,
            int[] retainedRows,
            int originalRows,
            ComputeBackend backend,
            BackendProvenance provenance) {
        boolean allowMinimumNorm = options.rankDeficiencyStrategy()
            == RankDeficiencyStrategy.MINIMUM_NORM;
        double[] coefficients = initialCoefficients(
            response, design, rows, columns, family,
            priorWeights, offset, options, allowMinimumNorm, backend);
        State state = state(response, design, rows, columns,
            family, priorWeights, offset, coefficients, backend);
        boolean converged = false;
        String message = "maximum iterations reached";
        int iterations = 0;

        for (int iteration = 1; iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            WorkingData working = workingData(
                response, design, rows, columns, family,
                priorWeights, offset, state);
            Solution solution = LeastSquaresSolver.solve(
                working.design(), working.response(), rows, columns,
                allowMinimumNorm, backend);
            double[] candidateCoefficients = solution.coefficients();
            State candidate = null;
            double scale = 1.0;
            for (int attempt = 0; attempt < 30; attempt++) {
                double[] trial = interpolate(
                    coefficients, candidateCoefficients, scale);
                State trialState = state(response, design, rows, columns,
                    family, priorWeights, offset, trial, backend);
                if (trialState.deviance() <= state.deviance()
                        + 1e-12 * (1.0 + state.deviance())) {
                    candidateCoefficients = trial;
                    candidate = trialState;
                    break;
                }
                scale *= 0.5;
            }
            if (candidate == null) {
                message = "IRLS step-halving could not reduce deviance";
                break;
            }

            double coefficientChange = relativeMaximumChange(
                coefficients, candidateCoefficients);
            double devianceChange = Math.abs(candidate.deviance() - state.deviance())
                / (1.0 + Math.abs(state.deviance()));
            coefficients = candidateCoefficients;
            state = candidate;
            if (coefficientChange <= options.relativeTolerance()
                    && devianceChange <= options.relativeTolerance()) {
                converged = true;
                message = "coefficient and deviance tolerances reached";
                break;
            }
        }

        WorkingData finalWorking = workingData(
            response, design, rows, columns, family,
            priorWeights, offset, state);
        Solution finalSolution = LeastSquaresSolver.solve(
            finalWorking.design(), finalWorking.response(), rows, columns,
            allowMinimumNorm, backend);
        int degreesOfFreedom = rows - finalSolution.rank();
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException(
                "GLM requires at least one residual degree of freedom");
        }

        double pearsonChiSquare = 0.0;
        double[] devianceResiduals = new double[rows];
        double[] pearsonResiduals = new double[rows];
        for (int row = 0; row < rows; row++) {
            double mean = state.means()[row];
            double variance = family.variance(mean);
            double residual = response[row] - mean;
            pearsonResiduals[row] = residual
                * Math.sqrt(priorWeights[row] / variance);
            pearsonChiSquare += pearsonResiduals[row] * pearsonResiduals[row];
            devianceResiduals[row] = Math.copySign(
                Math.sqrt(Math.max(0.0, priorWeights[row]
                    * family.unitDeviance(response[row], mean))), residual);
        }

        boolean estimateDispersion = !family.fixedDispersion()
            || options.dispersionMode() == DispersionMode.PEARSON;
        double dispersion = estimateDispersion
            ? pearsonChiSquare / degreesOfFreedom : 1.0;
        dispersion = Math.max(Double.MIN_NORMAL, dispersion);
        double[] covariance = finalSolution.unscaledCovariance().clone();
        for (int index = 0; index < covariance.length; index++) {
            covariance[index] *= dispersion;
        }

        double[] standardErrors = new double[columns];
        double[] statistics = new double[columns];
        double[] pValues = new double[columns];
        double[] confidenceLower = new double[columns];
        double[] confidenceUpper = new double[columns];
        double alpha = 1.0 - options.confidenceLevel();
        double critical = Normal.quantile(
            1.0 - alpha / 2.0, 0.0, 1.0, true, false);
        for (int column = 0; column < columns; column++) {
            standardErrors[column] = Math.sqrt(Math.max(0.0,
                covariance[column * columns + column]));
            if (standardErrors[column] == 0.0) {
                statistics[column] = coefficients[column] == 0.0
                    ? Double.NaN
                    : Math.copySign(Double.POSITIVE_INFINITY, coefficients[column]);
            } else {
                statistics[column] = coefficients[column] / standardErrors[column];
            }
            pValues[column] = Double.isNaN(statistics[column])
                ? Double.NaN
                : Math.min(1.0, 2.0 * Normal.cumulative(
                    Math.abs(statistics[column]), 0.0, 1.0, false, false));
            double margin = critical * standardErrors[column];
            confidenceLower[column] = coefficients[column] - margin;
            confidenceUpper[column] = coefficients[column] + margin;
        }

        double logLikelihood = 0.0;
        for (int row = 0; row < rows; row++) {
            double contribution = family.logLikelihood(response[row],
                state.means()[row], priorWeights[row], dispersion);
            if (Double.isNaN(contribution)) {
                logLikelihood = Double.NaN;
                break;
            }
            logLikelihood += contribution;
        }
        int likelihoodParameterCount = columns + (estimateDispersion ? 1 : 0);
        double aic = Double.isNaN(logLikelihood) ? Double.NaN
            : 2.0 * likelihoodParameterCount - 2.0 * logLikelihood;

        return new GlmResult(
            family.name(), coefficients, covariance,
            standardErrors, statistics, pValues,
            confidenceLower, confidenceUpper,
            state.linearPredictor(), state.means(),
            devianceResiduals, pearsonResiduals,
            state.deviance(), dispersion, logLikelihood, aic,
            rows, columns, finalSolution.rank(), degreesOfFreedom,
            iterations, converged, message,
            retainedRows, originalRows, provenance);
    }

    private static double[] initialCoefficients(
            double[] response, double[] design, int rows, int columns,
            GlmFamily family, double[] priorWeights, double[] offset,
            GlmOptions options, boolean allowMinimumNorm,
            ComputeBackend backend) {
        double[] supplied = options.initialCoefficients();
        if (supplied != null) {
            if (supplied.length != columns) {
                throw new IllegalArgumentException(
                    "initial coefficient count must equal design columns");
            }
            return supplied;
        }
        double[] target = new double[rows];
        for (int row = 0; row < rows; row++) {
            target[row] = family.link(family.initialMean(response[row])) - offset[row];
        }
        WorkingData weighted = weight(design, target, rows, columns, priorWeights);
        return LeastSquaresSolver.solve(weighted.design(), weighted.response(),
            rows, columns, allowMinimumNorm, backend).coefficients();
    }

    private static State state(
            double[] response, double[] design, int rows, int columns,
            GlmFamily family, double[] priorWeights, double[] offset,
            double[] coefficients, ComputeBackend backend) {
        double[] predictor = MatrixOps.multiply(
            backend, design, rows, columns, coefficients);
        double[] means = new double[rows];
        double deviance = 0.0;
        for (int row = 0; row < rows; row++) {
            predictor[row] += offset[row];
            means[row] = family.inverseLink(predictor[row]);
            if (!Double.isFinite(means[row])) {
                throw new IllegalArgumentException(
                    "inverse link produced a non-finite mean");
            }
            deviance += priorWeights[row]
                * family.unitDeviance(response[row], means[row]);
        }
        return new State(predictor, means, deviance);
    }

    private static WorkingData workingData(
            double[] response, double[] design, int rows, int columns,
            GlmFamily family, double[] priorWeights, double[] offset,
            State state) {
        double[] targets = new double[rows];
        double[] weights = new double[rows];
        for (int row = 0; row < rows; row++) {
            double derivative = family.meanDerivative(state.linearPredictor()[row]);
            double variance = family.variance(state.means()[row]);
            if (!Double.isFinite(derivative) || derivative == 0.0
                    || !Double.isFinite(variance) || variance <= 0.0) {
                throw new IllegalArgumentException(
                    "family produced invalid IRLS derivative or variance");
            }
            weights[row] = clamp(priorWeights[row] * derivative * derivative / variance,
                MINIMUM_WORKING_WEIGHT, MAXIMUM_WORKING_WEIGHT);
            targets[row] = state.linearPredictor()[row]
                + (response[row] - state.means()[row]) / derivative
                - offset[row];
        }
        return weight(design, targets, rows, columns, weights);
    }

    private static WorkingData weight(
            double[] design, double[] response, int rows, int columns,
            double[] weights) {
        double[] weightedDesign = new double[design.length];
        double[] weightedResponse = new double[rows];
        for (int row = 0; row < rows; row++) {
            double scale = Math.sqrt(weights[row]);
            weightedResponse[row] = scale * response[row];
            for (int column = 0; column < columns; column++) {
                weightedDesign[row * columns + column] =
                    scale * design[row * columns + column];
            }
        }
        return new WorkingData(weightedDesign, weightedResponse);
    }

    private static double[] prepareWeights(double[] weights, int rows) {
        if (weights == null) {
            double[] result = new double[rows];
            Arrays.fill(result, 1.0);
            return result;
        }
        if (weights.length != rows) {
            throw new IllegalArgumentException("prior weight length must equal rows");
        }
        double[] result = MatrixOps.finiteCopy(weights, "priorWeights");
        for (double value : result) {
            if (!(value > 0.0)) {
                throw new IllegalArgumentException(
                    "prior weights must be strictly positive");
            }
        }
        return result;
    }

    private static double[] prepareOffset(double[] offset, int rows) {
        if (offset == null) {
            return new double[rows];
        }
        if (offset.length != rows) {
            throw new IllegalArgumentException("offset length must equal rows");
        }
        return MatrixOps.finiteCopy(offset, "offset");
    }

    private static void validateResponses(
            double[] response, double[] weights, GlmFamily family) {
        for (int row = 0; row < response.length; row++) {
            family.validateResponse(response[row], weights[row]);
        }
    }

    private static double[] interpolate(
            double[] current, double[] candidate, double scale) {
        double[] result = new double[current.length];
        for (int index = 0; index < current.length; index++) {
            result[index] = current[index]
                + scale * (candidate[index] - current[index]);
        }
        return result;
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

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    private record State(
            double[] linearPredictor,
            double[] means,
            double deviance) {
    }

    private record WorkingData(double[] design, double[] response) {
    }
}
