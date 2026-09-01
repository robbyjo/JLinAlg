/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

import java.util.ArrayList;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;

/** Gaussian REML for dense positive-semidefinite covariance components. */
public final class Reml {
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

    private Reml() {
    }

    /** Fits REML using default options and the preferred compute policy. */
    public static RemlResult fit(
            double[] response,
            double[][] fixedEffects,
            List<VarianceComponent> components) {
        return fit(response, fixedEffects, components,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits REML from a conventional rectangular fixed-effect matrix. */
    public static RemlResult fit(
            double[] response,
            double[][] fixedEffects,
            List<VarianceComponent> components,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, rowMajor, response.length, fixedEffects[0].length,
            components, options, backendPolicy);
    }

    /** Fits REML from a contiguous row-major fixed-effect matrix. */
    public static RemlResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<VarianceComponent> components,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        validateComponents(components, rows);
        if (rows <= columns) {
            throw new IllegalArgumentException(
                "REML requires more observations than fixed-effect columns");
        }
        if (options == null) {
            throw new IllegalArgumentException("options are required");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fit(response, fixedEffects, rows, columns, components,
                new double[rows * rows], options,
                context.backend(), context.provenance());
        }
    }

    /**
     * Fits REML with a known, unscaled covariance contribution in addition to
     * the covariance components whose positive scales are estimated.
     */
    public static RemlResult fitWithKnownCovariance(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<VarianceComponent> components,
            double[] knownCovariance,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (backendPolicy == null) {
            throw new IllegalArgumentException("backendPolicy is required");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            return fitWithKnownCovariance(response, fixedEffects, rows, columns,
                components, knownCovariance, options, context);
        }
    }

    /** Fits REML with known covariance using a caller-owned backend context. */
    public static RemlResult fitWithKnownCovariance(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<VarianceComponent> components,
            double[] knownCovariance,
            RemlOptions options,
            BackendContext context) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        validateComponents(components, rows);
        if (rows <= columns) {
            throw new IllegalArgumentException(
                "REML requires more observations than fixed-effect columns");
        }
        if (options == null || context == null) {
            throw new IllegalArgumentException(
                "options and context are required");
        }
        double[] known = new VarianceComponent(
            "known covariance", rows, knownCovariance).covarianceView().clone();
        return fit(response, fixedEffects, rows, columns, components,
            known, options, context.backend(), context.provenance());
    }

    private static RemlResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<VarianceComponent> components,
            double[] knownCovariance,
            RemlOptions options,
            ComputeBackend backend,
            BackendProvenance provenance) {
        double[] initial = initialVariances(response, components, options);
        double minimumLog = Math.log(options.minimumVariance());
        double maximumLog = Math.log(options.maximumVariance());
        double[] logVariances = new double[initial.length];
        for (int index = 0; index < initial.length; index++) {
            logVariances[index] = clamp(
                Math.log(initial[index]), minimumLog, maximumLog);
        }

        Evaluation current = evaluate(response, fixedEffects, rows, columns,
            components, knownCovariance, logVariances,
            options.varianceEstimation(), backend);
        boolean converged = false;
        String message = "maximum iterations reached";
        int iterations = 0;

        for (int iteration = 1; iteration <= options.maximumIterations(); iteration++) {
            iterations = iteration;
            double projectedScore = projectedScoreMaximum(
                current.score(), logVariances, minimumLog, maximumLog);
            if (projectedScore <= options.scoreTolerance()) {
                converged = true;
                message = "projected score tolerance reached";
                break;
            }

            double[] step = solveInformation(
                current.information(), current.score(), backend);
            limitMaximumAbsolute(step, options.maximumLogVarianceStep());
            SearchResult search = lineSearch(
                response, fixedEffects, rows, columns, components,
                knownCovariance, logVariances, current, step,
                minimumLog, maximumLog, options.varianceEstimation(), backend);

            if (search == null) {
                double[] gradientStep = current.score().clone();
                limitMaximumAbsolute(
                    gradientStep, options.maximumLogVarianceStep());
                search = lineSearch(
                    response, fixedEffects, rows, columns, components,
                    knownCovariance, logVariances, current, gradientStep,
                    minimumLog, maximumLog, options.varianceEstimation(), backend);
            }

            if (search == null) {
                message = "line search could not improve the likelihood";
                break;
            }

            double relativeChange = Math.abs(
                search.evaluation().logLikelihood() - current.logLikelihood())
                / (1.0 + Math.abs(current.logLikelihood()));
            logVariances = search.logVariances();
            current = search.evaluation();

            if (relativeChange <= options.relativeTolerance()
                    && projectedScoreMaximum(current.score(), logVariances,
                        minimumLog, maximumLog) <= options.scoreTolerance()) {
                converged = true;
                message = "likelihood and projected score tolerances reached";
                break;
            }
        }

        double[] fitted = MatrixOps.multiply(
            backend, fixedEffects, rows, columns, current.fixedEffects());
        double[] residuals = MatrixOps.subtract(response, fitted);
        List<String> names = new ArrayList<>(components.size());
        for (VarianceComponent component : components) {
            names.add(component.name());
        }
        FixedEffectInference inference = fixedEffectInference(
            current, components, rows, columns,
            options.degreesOfFreedomMethod(), backend);
        double[] standardErrors = new double[columns];
        for (int column = 0; column < columns; column++) {
            standardErrors[column] = Math.sqrt(Math.max(0.0,
                inference.covariance()[column * columns + column]));
        }
        AssociationStatistics association = AssociationStatistics.studentT(
            current.fixedEffects(), standardErrors, inference.degreesOfFreedom(),
            options.degreesOfFreedomMethod());

        return new RemlResult(
            names, current.variances(), current.fixedEffects(),
            current.fixedEffectCovariance(), inference.covariance(), standardErrors,
            fitted, residuals, current.score(), association,
            current.logLikelihood(),
            options.varianceEstimation(),
            rows, columns, iterations, converged, message, provenance);
    }

    private static Evaluation evaluate(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<VarianceComponent> components,
            double[] knownCovariance,
            double[] logVariances,
            VarianceEstimation estimation,
            ComputeBackend backend) {
        int componentCount = components.size();
        double[] variances = new double[componentCount];
        double[] covariance = knownCovariance.clone();
        for (int component = 0; component < componentCount; component++) {
            variances[component] = Math.exp(logVariances[component]);
            double[] basis = components.get(component).covarianceView();
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] += variances[component] * basis[index];
            }
        }

        CholeskyFactor covarianceFactor;
        CholeskyFactor fixedFactor;
        try {
            covarianceFactor = backend.dpotrf(covariance, rows);
            double[] inverseCovarianceFixed =
                covarianceFactor.solve(fixedEffects, columns);
            double[] fixedInformation = MatrixOps.transposeMultiply(
                backend, fixedEffects, rows, columns,
                inverseCovarianceFixed, columns);
            symmetrize(fixedInformation, columns);
            fixedFactor = backend.dpotrf(fixedInformation, columns);

            double[] inverseCovarianceResponse = covarianceFactor.solve(response);
            double[] fixedRightSide = new double[columns];
            backend.dgemv(MatrixTranspose.TRANSPOSE, rows, columns,
                1.0, fixedEffects, inverseCovarianceResponse,
                0.0, fixedRightSide);
            double[] beta = fixedFactor.solve(fixedRightSide);
            double[] residual = MatrixOps.subtract(response,
                MatrixOps.multiply(backend, fixedEffects, rows, columns, beta));
            double[] projectedResponse = covarianceFactor.solve(residual);
            double quadratic = backend.ddot(rows, residual, 0, 1,
                projectedResponse, 0, 1);
            boolean restricted = estimation == VarianceEstimation.REML;
            double logLikelihood = -0.5 * (
                (restricted ? rows - columns : rows) * LOG_TWO_PI
                + covarianceFactor.logDeterminant()
                + (restricted ? fixedFactor.logDeterminant() : 0.0)
                + quadratic);

            double[] fixedCovariance = fixedFactor.solve(
                MatrixOps.identity(columns), columns);
            symmetrize(fixedCovariance, columns);
            double[] inverseCovariance = covarianceFactor.solve(
                MatrixOps.identity(rows), rows);
            double[] temporary = MatrixOps.multiply(
                backend, inverseCovarianceFixed, rows, columns,
                fixedCovariance, columns);
            double[] correction = new double[rows * rows];
            backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                rows, rows, columns, 1.0,
                temporary, inverseCovarianceFixed, 0.0, correction);
            double[] projection = MatrixOps.subtract(
                inverseCovariance, correction);
            symmetrize(projection, rows);

            double[][] derivativeProducts = new double[componentCount][];
            double[] score = new double[componentCount];
            for (int component = 0; component < componentCount; component++) {
                double[] basis = components.get(component).covarianceView();
                derivativeProducts[component] = MatrixOps.multiply(
                    backend, restricted ? projection : inverseCovariance,
                    rows, rows, basis, rows);
                double trace = trace(derivativeProducts[component], rows);
                double[] basisTimesProjected = MatrixOps.multiply(
                    backend, basis, rows, rows, projectedResponse);
                double quadraticDerivative = backend.ddot(rows,
                    projectedResponse, 0, 1,
                    basisTimesProjected, 0, 1);
                score[component] = 0.5 * variances[component]
                    * (quadraticDerivative - trace);
            }

            double[] information = new double[componentCount * componentCount];
            for (int row = 0; row < componentCount; row++) {
                for (int column = 0; column <= row; column++) {
                    double value = 0.5 * variances[row] * variances[column]
                        * traceProduct(derivativeProducts[row],
                            derivativeProducts[column], rows);
                    information[row * componentCount + column] = value;
                    information[column * componentCount + row] = value;
                }
            }
            return new Evaluation(variances, beta, fixedCovariance,
                inverseCovariance, inverseCovarianceFixed, temporary,
                score, information, logLikelihood);
        } catch (IllegalArgumentException | IllegalStateException exception) {
            throw new IllegalArgumentException(
                "covariance or fixed-effect information matrix is not positive definite",
                exception);
        }
    }

    private static SearchResult lineSearch(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<VarianceComponent> components,
            double[] knownCovariance,
            double[] currentLogVariances,
            Evaluation current,
            double[] direction,
            double minimumLog,
            double maximumLog,
            VarianceEstimation estimation,
            ComputeBackend backend) {
        double scale = 1.0;
        for (int attempt = 0; attempt < 24; attempt++) {
            double[] candidate = currentLogVariances.clone();
            double maximumChange = 0.0;
            for (int component = 0; component < candidate.length; component++) {
                double updated = clamp(candidate[component]
                    + scale * direction[component], minimumLog, maximumLog);
                maximumChange = Math.max(maximumChange,
                    Math.abs(updated - candidate[component]));
                candidate[component] = updated;
            }
            if (maximumChange <= 1e-14) {
                return null;
            }
            try {
                Evaluation evaluated = evaluate(response, fixedEffects,
                    rows, columns, components, knownCovariance, candidate,
                    estimation, backend);
                if (evaluated.logLikelihood() > current.logLikelihood()) {
                    return new SearchResult(candidate, evaluated);
                }
            } catch (IllegalArgumentException ignored) {
                // A smaller positive-definite step may still be valid.
            }
            scale *= 0.5;
        }
        return null;
    }

    private static double[] solveInformation(
            double[] information, double[] score, ComputeBackend backend) {
        int dimension = score.length;
        double maximumDiagonal = 0.0;
        for (int index = 0; index < dimension; index++) {
            maximumDiagonal = Math.max(maximumDiagonal,
                Math.abs(information[index * dimension + index]));
        }
        double ridge = Math.max(1e-12, maximumDiagonal * 1e-10);
        for (int attempt = 0; attempt < 10; attempt++) {
            double[] regularized = information.clone();
            for (int index = 0; index < dimension; index++) {
                regularized[index * dimension + index] += ridge;
            }
            try {
                return backend.dpotrf(regularized, dimension).solve(score);
            } catch (IllegalArgumentException | IllegalStateException ignored) {
                ridge *= 10.0;
            }
        }
        throw new IllegalArgumentException(
            "variance components are not numerically identifiable");
    }

    private static FixedEffectInference fixedEffectInference(
            Evaluation evaluation,
            List<VarianceComponent> components,
            int rows,
            int columns,
            DegreesOfFreedomMethod method,
            ComputeBackend backend) {
        double[] result = new double[columns];
        if (method == DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION) {
            double degrees = rows - columns - 1.0;
            for (int column = 0; column < columns; column++) {
                result[column] = degrees;
            }
            return new FixedEffectInference(
                evaluation.fixedEffectCovariance(), result);
        }
        if (method == DegreesOfFreedomMethod.KENWARD_ROGER) {
            double[] logVarianceCovariance = inverseInformation(
                evaluation.information(), components.size(), backend,
                "Kenward-Roger adjustment");
            KenwardRoger.Result adjusted = KenwardRoger.calculate(
                evaluation.fixedEffectCovariance(),
                evaluation.inverseCovariance(),
                evaluation.inverseCovarianceFixed(),
                logVarianceCovariance, evaluation.variances(), components,
                rows, columns, backend);
            return new FixedEffectInference(
                adjusted.adjustedCovariance(), adjusted.degreesOfFreedom());
        }
        if (method != DegreesOfFreedomMethod.SATTERTHWAITE) {
            throw new IllegalArgumentException(
                "unsupported REML degrees-of-freedom method: " + method);
        }

        int componentCount = components.size();
        double[] varianceParameterCovariance = inverseInformation(
            evaluation.information(), componentCount, backend,
            "Satterthwaite DF");
        double[] gradients = new double[componentCount * columns];
        double[] covarianceProjection = evaluation.covarianceProjection();
        for (int component = 0; component < componentCount; component++) {
            double[] basisTimesProjection = MatrixOps.multiply(
                backend, components.get(component).covarianceView(),
                rows, rows, covarianceProjection, columns);
            double variance = evaluation.variances()[component];
            for (int coefficient = 0; coefficient < columns; coefficient++) {
                double derivative = 0.0;
                for (int row = 0; row < rows; row++) {
                    derivative += covarianceProjection[row * columns + coefficient]
                        * basisTimesProjection[row * columns + coefficient];
                }
                gradients[component * columns + coefficient] =
                    variance * derivative;
            }
        }

        double[] fixedCovariance = evaluation.fixedEffectCovariance();
        for (int coefficient = 0; coefficient < columns; coefficient++) {
            double variance = fixedCovariance[coefficient * columns + coefficient];
            double varianceOfVariance = 0.0;
            for (int left = 0; left < componentCount; left++) {
                double leftGradient = gradients[left * columns + coefficient];
                for (int right = 0; right < componentCount; right++) {
                    varianceOfVariance += leftGradient
                        * varianceParameterCovariance[left * componentCount + right]
                        * gradients[right * columns + coefficient];
                }
            }
            result[coefficient] = variance > 0.0
                    && varianceOfVariance > 0.0
                    && Double.isFinite(varianceOfVariance)
                ? 2.0 * variance * variance / varianceOfVariance
                : Double.NaN;
        }
        return new FixedEffectInference(
            evaluation.fixedEffectCovariance(), result);
    }

    private static double[] inverseInformation(
            double[] information,
            int dimension,
            ComputeBackend backend,
            String calculation) {
        double maximumDiagonal = 0.0;
        for (int index = 0; index < dimension; index++) {
            maximumDiagonal = Math.max(maximumDiagonal,
                Math.abs(information[index * dimension + index]));
        }
        double ridge = Math.max(1e-12, maximumDiagonal * 1e-10);
        for (int attempt = 0; attempt < 10; attempt++) {
            double[] regularized = information.clone();
            for (int index = 0; index < dimension; index++) {
                regularized[index * dimension + index] += ridge;
            }
            try {
                double[] inverse = backend.dpotrf(regularized, dimension).solve(
                    MatrixOps.identity(dimension), dimension);
                symmetrize(inverse, dimension);
                return inverse;
            } catch (IllegalArgumentException | IllegalStateException ignored) {
                ridge *= 10.0;
            }
        }
        throw new IllegalArgumentException(
            "variance-component information is not invertible for " + calculation);
    }

    private static double[] initialVariances(
            double[] response,
            List<VarianceComponent> components,
            RemlOptions options) {
        double[] supplied = options.initialVariances();
        if (supplied != null) {
            if (supplied.length != components.size()) {
                throw new IllegalArgumentException(
                    "initial variance count must equal component count");
            }
            for (double value : supplied) {
                if (value < options.minimumVariance()
                        || value > options.maximumVariance()) {
                    throw new IllegalArgumentException(
                        "initial variances must lie within the configured bounds");
                }
            }
            return supplied;
        }

        double mean = 0.0;
        for (double value : response) {
            mean += value;
        }
        mean /= response.length;
        double responseVariance = 0.0;
        for (double value : response) {
            double centered = value - mean;
            responseVariance += centered * centered;
        }
        responseVariance /= Math.max(1, response.length - 1);
        responseVariance = Math.max(responseVariance, options.minimumVariance());

        double[] result = new double[components.size()];
        for (int component = 0; component < components.size(); component++) {
            VarianceComponent value = components.get(component);
            double[] basis = value.covarianceView();
            double meanDiagonal = 0.0;
            for (int index = 0; index < value.dimension(); index++) {
                meanDiagonal += basis[index * value.dimension() + index];
            }
            meanDiagonal /= value.dimension();
            double estimate = responseVariance / components.size()
                / Math.max(meanDiagonal, 1e-12);
            result[component] = clamp(estimate,
                options.minimumVariance(), options.maximumVariance());
        }
        return result;
    }

    private static void validateComponents(
            List<VarianceComponent> components, int dimension) {
        if (components == null || components.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one variance component is required");
        }
        for (VarianceComponent component : components) {
            if (component == null || component.dimension() != dimension) {
                throw new IllegalArgumentException(
                    "all variance components must match the response dimension");
            }
        }
    }

    private static double projectedScoreMaximum(
            double[] score,
            double[] logVariances,
            double minimumLog,
            double maximumLog) {
        double maximum = 0.0;
        for (int index = 0; index < score.length; index++) {
            double value = score[index];
            if ((logVariances[index] <= minimumLog + 1e-12 && value < 0.0)
                    || (logVariances[index] >= maximumLog - 1e-12 && value > 0.0)) {
                value = 0.0;
            }
            maximum = Math.max(maximum, Math.abs(value));
        }
        return maximum;
    }

    private static void limitMaximumAbsolute(double[] values, double limit) {
        double maximum = 0.0;
        for (double value : values) {
            maximum = Math.max(maximum, Math.abs(value));
        }
        if (maximum > limit) {
            double scale = limit / maximum;
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

    private static double trace(double[] matrix, int dimension) {
        double result = 0.0;
        for (int index = 0; index < dimension; index++) {
            result += matrix[index * dimension + index];
        }
        return result;
    }

    private static double traceProduct(
            double[] left, double[] right, int dimension) {
        double result = 0.0;
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                result += left[row * dimension + column]
                    * right[column * dimension + row];
            }
        }
        return result;
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    private record Evaluation(
            double[] variances,
            double[] fixedEffects,
            double[] fixedEffectCovariance,
            double[] inverseCovariance,
            double[] inverseCovarianceFixed,
            double[] covarianceProjection,
            double[] score,
            double[] information,
            double logLikelihood) {
    }

    private record SearchResult(
            double[] logVariances,
            Evaluation evaluation) {
    }

    private record FixedEffectInference(
            double[] covariance,
            double[] degreesOfFreedom) {
    }
}
