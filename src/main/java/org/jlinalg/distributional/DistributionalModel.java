/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.Normal;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;
import org.jlinalg.internal.MatrixOps;

/** Multi-predictor penalized maximum likelihood using block Fisher scoring. */
public final class DistributionalModel {
    private DistributionalModel() { }

    /** Fits with default controls and preferred acceleration. */
    public static DistributionalResult fit(
            double[] response,
            List<PenalizedPredictor> predictors,
            DistributionalFamily family) {
        return fit(response, predictors, family,
            DistributionalOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a distributional or vector additive model. */
    public static DistributionalResult fit(
            double[] response,
            List<PenalizedPredictor> predictors,
            DistributionalFamily family,
            DistributionalOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, predictors, family, options, backendPolicy, null);
    }

    /** Fits with optional parameter-block coefficients used as warm starts. */
    public static DistributionalResult fit(
            double[] response,
            List<PenalizedPredictor> predictors,
            DistributionalFamily family,
            DistributionalOptions options,
            BackendPolicy backendPolicy,
            List<double[]> startingCoefficients) {
        if (response == null || response.length == 0
                || predictors == null || family == null
                || options == null || backendPolicy == null
                || predictors.size() != family.parameterCount()) {
            throw new IllegalArgumentException(
                "response, one predictor per family parameter, controls, and backend are required");
        }
        MatrixOps.requireFinite(response, "response");
        for (double value : response) family.validateResponse(value);
        int observations = response.length;
        int parameterCount = family.parameterCount();
        int[] starts = new int[parameterCount + 1];
        double[][] designs = new double[parameterCount][];
        double[][] penalties = new double[parameterCount][];
        for (int parameter = 0; parameter < parameterCount; parameter++) {
            PenalizedPredictor predictor = predictors.get(parameter);
            if (predictor == null || predictor.observations() != observations) {
                throw new IllegalArgumentException(
                    "all predictors must match the response length");
            }
            starts[parameter + 1] = starts[parameter] + predictor.columns();
            designs[parameter] = predictor.design();
            penalties[parameter] = predictor.penaltyDiagonal();
        }
        int totalColumns = starts[parameterCount];
        double[] coefficients = startingCoefficients == null
            ? initialize(response, predictors, family, starts, designs)
            : starting(startingCoefficients, predictors, starts);
        State current = evaluate(
            response, coefficients, starts, designs, family);
        boolean converged = false;
        String message = "maximum Fisher-scoring iterations reached";
        int iterations = 0;
        SystemState system = null;

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            for (int iteration = 1;
                    iteration <= options.maximumIterations(); iteration++) {
                iterations = iteration;
                system = system(response, coefficients, current, starts,
                    designs, penalties, family, totalColumns);
                double[] step;
                try {
                    step = backend.dpotrf(system.information(), totalColumns)
                        .solve(system.gradient());
                } catch (IllegalArgumentException exception) {
                    message = "Fisher information is not positive definite";
                    break;
                }
                limit(step, options.maximumStep());
                State candidate = null;
                double[] candidateCoefficients = null;
                double scale = 1.0;
                for (int attempt = 0; attempt < 30; attempt++) {
                    double[] trial = coefficients.clone();
                    for (int column = 0; column < totalColumns; column++) {
                        trial[column] += scale * step[column];
                    }
                    State trialState = evaluate(
                        response, trial, starts, designs, family);
                    double trialPenalized = trialState.logLikelihood()
                        - 0.5 * penaltyQuadratic(trial, starts, penalties);
                    if (Double.isFinite(trialPenalized)
                            && trialPenalized >= system.penalizedLogLikelihood()
                                - 1e-12 * (1.0
                                    + Math.abs(system.penalizedLogLikelihood()))) {
                        candidate = trialState;
                        candidateCoefficients = trial;
                        break;
                    }
                    scale *= 0.5;
                }
                if (candidate == null) {
                    message = "step halving could not improve penalized likelihood";
                    break;
                }
                double coefficientChange = relativeMaximumChange(
                    coefficients, candidateCoefficients);
                double objectiveChange = Math.abs(
                    candidate.logLikelihood() - current.logLikelihood())
                    / (1.0 + Math.abs(current.logLikelihood()));
                coefficients = candidateCoefficients;
                current = candidate;
                if (coefficientChange <= options.relativeTolerance()
                        && objectiveChange <= options.relativeTolerance()) {
                    converged = true;
                    message = "coefficient and likelihood tolerances reached";
                    break;
                }
            }

            system = system(response, coefficients, current, starts,
                designs, penalties, family, totalColumns);
            CholeskyFactor factor = backend.dpotrf(
                system.information(), totalColumns);
            double[] covariance = factor.solve(
                MatrixOps.identity(totalColumns), totalColumns);
            List<DistributionalParameterResult> parameterResults =
                parameterResults(coefficients, covariance, current,
                    predictors, starts, totalColumns, family);
            return new DistributionalResult(family.name(), parameterResults,
                covariance, current.logLikelihood(),
                system.penalizedLogLikelihood(), iterations,
                converged, message);
        }
    }

    private static double[] starting(
            List<double[]> supplied,
            List<PenalizedPredictor> predictors,
            int[] starts) {
        if (supplied.size() != predictors.size()) {
            throw new IllegalArgumentException(
                "one starting coefficient block is required per parameter");
        }
        double[] result = new double[starts[starts.length - 1]];
        for (int parameter = 0; parameter < supplied.size(); parameter++) {
            double[] values = supplied.get(parameter);
            if (values == null
                    || values.length != predictors.get(parameter).columns()) {
                throw new IllegalArgumentException(
                    "starting coefficient dimensions do not match predictors");
            }
            MatrixOps.requireFinite(values, "starting coefficients");
            System.arraycopy(values, 0, result, starts[parameter], values.length);
        }
        return result;
    }

    private static double[] initialize(
            double[] response,
            List<PenalizedPredictor> predictors,
            DistributionalFamily family,
            int[] starts,
            double[][] designs) {
        double[] initial = family.initialParameters(response);
        if (initial.length != family.parameterCount()) {
            throw new IllegalArgumentException(
                "family returned the wrong number of initial parameters");
        }
        double[] result = new double[starts[starts.length - 1]];
        for (int parameter = 0; parameter < initial.length; parameter++) {
            double target = family.link(parameter, initial[parameter]);
            int columns = predictors.get(parameter).columns();
            if (columns > 0) {
                double intercept = designs[parameter][0];
                boolean constant = intercept != 0.0;
                for (int row = 1; row < response.length; row++) {
                    constant &= designs[parameter][row * columns] == intercept;
                }
                if (constant) result[starts[parameter]] = target / intercept;
            }
        }
        return result;
    }

    private static State evaluate(
            double[] response,
            double[] coefficients,
            int[] starts,
            double[][] designs,
            DistributionalFamily family) {
        int parameters = family.parameterCount();
        int observations = response.length;
        double[][] fitted = new double[parameters][observations];
        double[] values = new double[parameters];
        double logLikelihood = 0.0;
        for (int row = 0; row < observations; row++) {
            for (int parameter = 0; parameter < parameters; parameter++) {
                int columns = starts[parameter + 1] - starts[parameter];
                double predictor = 0.0;
                for (int column = 0; column < columns; column++) {
                    predictor += designs[parameter][row * columns + column]
                        * coefficients[starts[parameter] + column];
                }
                values[parameter] = family.inverseLink(parameter, predictor);
                fitted[parameter][row] = values[parameter];
            }
            logLikelihood += family.logLikelihood(response[row], values);
        }
        return new State(fitted, logLikelihood);
    }

    private static SystemState system(
            double[] response,
            double[] coefficients,
            State state,
            int[] starts,
            double[][] designs,
            double[][] penalties,
            DistributionalFamily family,
            int totalColumns) {
        int parameters = family.parameterCount();
        double[] gradient = new double[totalColumns];
        double[] information = new double[totalColumns * totalColumns];
        double[] values = new double[parameters];
        double[] score = new double[parameters];
        double[] observationInformation = new double[parameters * parameters];
        for (int row = 0; row < response.length; row++) {
            for (int parameter = 0; parameter < parameters; parameter++) {
                values[parameter] = state.fittedParameters()[parameter][row];
            }
            family.derivatives(response[row], values,
                score, observationInformation);
            for (int first = 0; first < parameters; first++) {
                int firstColumns = starts[first + 1] - starts[first];
                for (int left = 0; left < firstColumns; left++) {
                    int leftIndex = starts[first] + left;
                    double leftValue = designs[first][row * firstColumns + left];
                    gradient[leftIndex] += leftValue * score[first];
                    for (int second = 0; second < parameters; second++) {
                        int secondColumns = starts[second + 1] - starts[second];
                        double weight = observationInformation[
                            first * parameters + second];
                        for (int right = 0; right < secondColumns; right++) {
                            int rightIndex = starts[second] + right;
                            information[leftIndex * totalColumns + rightIndex] +=
                                leftValue * weight
                                    * designs[second][row * secondColumns + right];
                        }
                    }
                }
            }
        }
        for (int parameter = 0; parameter < parameters; parameter++) {
            for (int column = 0; column < penalties[parameter].length; column++) {
                int index = starts[parameter] + column;
                double penalty = penalties[parameter][column];
                gradient[index] -= penalty * coefficients[index];
                information[index * totalColumns + index] += penalty;
            }
        }
        symmetrize(information, totalColumns);
        double penalized = state.logLikelihood()
            - 0.5 * penaltyQuadratic(coefficients, starts, penalties);
        return new SystemState(
            gradient, information, penalized);
    }

    private static List<DistributionalParameterResult> parameterResults(
            double[] coefficients,
            double[] covariance,
            State state,
            List<PenalizedPredictor> predictors,
            int[] starts,
            int totalColumns,
            DistributionalFamily family) {
        List<DistributionalParameterResult> result = new ArrayList<>();
        for (int parameter = 0; parameter < family.parameterCount(); parameter++) {
            int columns = starts[parameter + 1] - starts[parameter];
            double[] beta = Arrays.copyOfRange(
                coefficients, starts[parameter], starts[parameter + 1]);
            double[] standardErrors = new double[columns];
            double[] statistics = new double[columns];
            double[] pValues = new double[columns];
            for (int column = 0; column < columns; column++) {
                int index = starts[parameter] + column;
                standardErrors[column] = Math.sqrt(Math.max(0.0,
                    covariance[index * totalColumns + index]));
                statistics[column] = standardErrors[column] == 0.0
                    ? Math.copySign(Double.POSITIVE_INFINITY, beta[column])
                    : beta[column] / standardErrors[column];
                pValues[column] = Math.min(1.0, 2.0 * Normal.cumulative(
                    Math.abs(statistics[column]), 0.0, 1.0, false, false));
            }
            double penaltyTrace = 0.0;
            double[] penalty = predictors.get(parameter).penaltyDiagonal();
            for (int column = 0; column < columns; column++) {
                int index = starts[parameter] + column;
                penaltyTrace += covariance[index * totalColumns + index]
                    * penalty[column];
            }
            double edf = Math.max(0.0,
                Math.min(columns, columns - penaltyTrace));
            result.add(new DistributionalParameterResult(
                family.parameterNames().get(parameter), beta,
                standardErrors, statistics, pValues,
                state.fittedParameters()[parameter], edf));
        }
        return List.copyOf(result);
    }

    private static double penaltyQuadratic(
            double[] coefficients, int[] starts, double[][] penalties) {
        double result = 0.0;
        for (int parameter = 0; parameter < penalties.length; parameter++) {
            for (int column = 0; column < penalties[parameter].length; column++) {
                double coefficient = coefficients[starts[parameter] + column];
                result += penalties[parameter][column]
                    * coefficient * coefficient;
            }
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

    private static void limit(double[] values, double maximum) {
        double observed = 0.0;
        for (double value : values) observed = Math.max(observed, Math.abs(value));
        if (observed > maximum) {
            double scale = maximum / observed;
            for (int index = 0; index < values.length; index++) {
                values[index] *= scale;
            }
        }
    }

    private static void symmetrize(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++) {
            for (int column = row + 1; column < dimension; column++) {
                double value = 0.5 * (matrix[row * dimension + column]
                    + matrix[column * dimension + row]);
                matrix[row * dimension + column] = value;
                matrix[column * dimension + row] = value;
            }
        }
    }

    private record State(
            double[][] fittedParameters,
            double logLikelihood) { }

    private record SystemState(
            double[] gradient,
            double[] information,
            double penalizedLogLikelihood) { }
}
