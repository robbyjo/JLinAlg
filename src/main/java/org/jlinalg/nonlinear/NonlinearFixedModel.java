/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

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
            LeastSquaresSolver.Solution finalLinear;
            for (; iteration < options.maximumIterations(); iteration++) {
                Evaluation evaluation = evaluate(response, parameters, mean);
                finalLinear = LeastSquaresSolver.solve(evaluation.gradient,
                    evaluation.residual, response.length, parameters.length,
                    false, backend);
                double[] step = finalLinear.coefficients();
                if (stationary(evaluation, step, response, options)) {
                    converged = true; iteration++; break;
                }
                double acceptedStep = options.initialStep();
                double[] candidate = parameters.clone();
                double candidateObjective = Double.POSITIVE_INFINITY;
                for (int halving = 0; halving <= options.maximumStepHalvings(); halving++) {
                    for (int j = 0; j < candidate.length; j++)
                        candidate[j] = parameters[j] + acceptedStep * step[j];
                    try { candidateObjective = objective(response, candidate, mean, null); }
                    catch (IllegalArgumentException invalidTrial) {
                        candidateObjective = Double.POSITIVE_INFINITY;
                    }
                    if (candidateObjective < objective) break;
                    acceptedStep *= 0.5;
                }
                if (!(candidateObjective < objective)) { iteration++; break; }
                parameters = candidate.clone(); objective = candidateObjective;
            }
            Evaluation finalEvaluation = evaluate(response, parameters, mean);
            // Inference must use the Jacobian at the returned parameters, even
            // when the iteration budget ended immediately after an accepted step.
            finalLinear = LeastSquaresSolver.solve(finalEvaluation.gradient,
                    finalEvaluation.residual, response.length, parameters.length,
                    false, backend);
            converged = converged || stationary(finalEvaluation, finalLinear.coefficients(),
                response, options);
            double variance = objective / (response.length - parameters.length);
            double[] covariance = finalLinear.unscaledCovariance().clone();
            for (int i = 0; i < covariance.length; i++) covariance[i] *= variance;
            double[] errors = new double[parameters.length];
            for (int i = 0; i < errors.length; i++)
                errors[i] = Math.sqrt(Math.max(0.0, covariance[i * errors.length + i]));
            return new NonlinearFitResult(AssociationStatistics.normal(parameters, errors),
                finalEvaluation.fitted, finalEvaluation.residual, objective, variance,
                iteration, converged, context.provenance());
        }
    }

    static Evaluation evaluate(double[] response, double[] parameters,
            NonlinearMeanFunction mean) {
        int rows = response.length, columns = parameters.length;
        double[] fitted = new double[rows], residual = new double[rows];
        double[] gradient = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            NonlinearMeanFunction.Evaluation value = mean.evaluate(parameters, row);
            double[] rowGradient = value.gradient();
            if (rowGradient.length != columns)
                throw new IllegalArgumentException("nonlinear gradient width must match parameters");
            fitted[row] = value.value(); residual[row] = response[row] - fitted[row];
            if (!Double.isFinite(residual[row]))
                throw new IllegalArgumentException("nonlinear residual is not finite");
            System.arraycopy(rowGradient, 0, gradient, row * columns, columns);
        }
        return new Evaluation(fitted, residual, gradient);
    }

    static double objective(double[] response, double[] parameters,
            NonlinearMeanFunction mean, double[] ignored) {
        Evaluation evaluation = evaluate(response, parameters, mean);
        double result = 0.0; for (double value : evaluation.residual) result += value * value;
        if (!Double.isFinite(result)) throw new IllegalArgumentException("nonlinear objective overflow");
        return result;
    }

    static void validate(double[] response, double[] initial,
            NonlinearMeanFunction mean, NonlinearModelOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || response.length < 2 || initial == null
                || initial.length == 0 || response.length <= initial.length || mean == null || options == null
                || backendPolicy == null) throw new IllegalArgumentException(
                    "response, parameters, mean, options, and backend are required");
        for (double value : response) if (!Double.isFinite(value))
            throw new IllegalArgumentException("response must be finite");
        for (double value : initial) if (!Double.isFinite(value))
            throw new IllegalArgumentException("initial parameters must be finite");
    }

    /** Scale-invariant Gauss-Newton stationarity in fitted-response space. */
    static boolean stationary(Evaluation evaluation, double[] step, double[] response, NonlinearModelOptions options) {
        // The projected correction controls both relative fitted-response
        // change and the fractional decrease predicted by the quadratic model.
        double tolerance=Math.min(options.parameterTolerance(),Math.sqrt(options.objectiveTolerance()));
        double projected = 0, residual = 0;
        for (int i = 0; i < response.length; i++) {
            double value = 0;
            for (int j = 0; j < step.length; j++) value += evaluation.gradient[i*step.length+j]*step[j];
            projected = Math.hypot(projected, value);
            residual = Math.hypot(residual, evaluation.residual[i]);
        }
        // An absolute response-scale floor is not a stationarity certificate:
        // an irrelevant large offset can otherwise hide a meaningful correction.
        // Precision-limited plateaus must not be reported as converged.
        return Double.isFinite(projected) && projected <= tolerance*residual;
    }

    record Evaluation(double[] fitted, double[] residual, double[] gradient) { }
}
