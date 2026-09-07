/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

import java.util.Arrays;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.LeastSquaresSolver;

/** Analytic-gradient Gauss-Newton nonlinear Gaussian fixed-effect fitting. */
public final class NonlinearFixedModel {
    private NonlinearFixedModel() { }

    public static NonlinearFitResult fit(double[] response,
            double[] initialParameters, NonlinearMeanFunction mean,
            NonlinearModelOptions options, BackendPolicy backendPolicy) {
        validate(response, initialParameters, mean, options, backendPolicy);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] parameters = initialParameters.clone();
            double objective = objective(response, parameters, mean, null);
            boolean converged = false; int iteration = 0;
            LeastSquaresSolver.Solution finalLinear = null;
            for (; iteration < options.maximumIterations(); iteration++) {
                Evaluation evaluation = evaluate(response, parameters, mean);
                finalLinear = LeastSquaresSolver.solve(evaluation.gradient,
                    evaluation.residual, response.length, parameters.length,
                    false, backend);
                double[] step = finalLinear.coefficients();
                double stepLength = norm(step);
                if (stepLength <= options.parameterTolerance()
                        * (1.0 + norm(parameters))) { converged = true; break; }
                double acceptedStep = options.initialStep();
                double[] candidate = parameters.clone();
                double candidateObjective = Double.POSITIVE_INFINITY;
                for (int halving = 0; halving <= options.maximumStepHalvings(); halving++) {
                    for (int j = 0; j < candidate.length; j++)
                        candidate[j] = parameters[j] + acceptedStep * step[j];
                    candidateObjective = objective(response, candidate, mean, null);
                    if (candidateObjective <= objective) break;
                    acceptedStep *= 0.5;
                }
                if (!(candidateObjective <= objective)) break;
                double improvement = objective - candidateObjective;
                parameters = candidate.clone(); objective = candidateObjective;
                if (improvement <= options.objectiveTolerance()
                        * (1.0 + objective)) { converged = true; break; }
            }
            Evaluation finalEvaluation = evaluate(response, parameters, mean);
            if (finalLinear == null)
                finalLinear = LeastSquaresSolver.solve(finalEvaluation.gradient,
                    finalEvaluation.residual, response.length, parameters.length,
                    false, backend);
            double variance = objective / Math.max(1, response.length - parameters.length);
            double[] covariance = finalLinear.unscaledCovariance().clone();
            for (int i = 0; i < covariance.length; i++) covariance[i] *= variance;
            double[] errors = new double[parameters.length];
            for (int i = 0; i < errors.length; i++)
                errors[i] = Math.sqrt(Math.max(0.0, covariance[i * errors.length + i]));
            return new NonlinearFitResult(AssociationStatistics.normal(parameters, errors),
                finalEvaluation.fitted, finalEvaluation.residual, objective, variance,
                iteration + 1, converged, context.provenance());
        }
    }

    static Evaluation evaluate(double[] response, double[] parameters,
            NonlinearMeanFunction mean) {
        int rows = response.length, columns = parameters.length;
        double[] fitted = new double[rows], residual = new double[rows];
        double[] gradient = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            NonlinearMeanFunction.Evaluation value = mean.evaluate(parameters, row);
            if (value.gradient().length != columns)
                throw new IllegalArgumentException("nonlinear gradient width must match parameters");
            fitted[row] = value.value(); residual[row] = response[row] - fitted[row];
            System.arraycopy(value.gradient(), 0, gradient, row * columns, columns);
        }
        return new Evaluation(fitted, residual, gradient);
    }

    static double objective(double[] response, double[] parameters,
            NonlinearMeanFunction mean, double[] ignored) {
        Evaluation evaluation = evaluate(response, parameters, mean);
        double result = 0.0; for (double value : evaluation.residual) result += value * value;
        return result;
    }

    private static void validate(double[] response, double[] initial,
            NonlinearMeanFunction mean, NonlinearModelOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || response.length < 2 || initial == null
                || initial.length == 0 || mean == null || options == null
                || backendPolicy == null) throw new IllegalArgumentException(
                    "response, parameters, mean, options, and backend are required");
        for (double value : response) if (!Double.isFinite(value))
            throw new IllegalArgumentException("response must be finite");
    }

    private static double norm(double[] values) {
        double result = 0.0; for (double value : values) result += value * value;
        return Math.sqrt(result);
    }

    record Evaluation(double[] fitted, double[] residual, double[] gradient) { }
}
