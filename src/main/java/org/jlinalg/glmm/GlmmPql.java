/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceComponent;

/**
 * Generalized linear mixed models fitted by first-order penalized
 * quasi-likelihood with REML working-model updates.
 */
public final class GlmmPql {
    private static final double MINIMUM_WORKING_WEIGHT = 1e-12;
    private static final double MAXIMUM_WORKING_WEIGHT = 1e150;

    private GlmmPql() { }

    /** Fits a PQL GLMM with unit prior weights, zero offset, and default controls. */
    public static GlmmPqlResult fit(
            double[] response,
            double[][] fixedEffects,
            GlmFamily family,
            List<VarianceComponent> randomComponents) {
        return fit(response, fixedEffects, family, randomComponents,
            null, null, GlmmPqlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits a PQL GLMM from a conventional fixed-effect matrix. */
    public static GlmmPqlResult fit(
            double[] response,
            double[][] fixedEffects,
            GlmFamily family,
            List<VarianceComponent> randomComponents,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, rowMajor, response.length, fixedEffects[0].length,
            family, randomComponents, priorWeights, offset, options, backendPolicy);
    }

    /** Fits a PQL GLMM from a contiguous row-major fixed-effect matrix. */
    public static GlmmPqlResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            GlmFamily family,
            List<VarianceComponent> randomComponents,
            double[] priorWeights,
            double[] offset,
            GlmmPqlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        if (family == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "family, options, and backendPolicy are required");
        }
        if (!family.fixedDispersion()) {
            throw new IllegalArgumentException(
                "PQL is provided for fixed-dispersion non-Gaussian families; "
                + "use Gaussian REML directly for gaussian(identity)");
        }
        validateComponents(randomComponents, rows);
        double[] weights = prepareWeights(priorWeights, rows);
        double[] offsets = prepareOffset(offset, rows);
        for (int row = 0; row < rows; row++) {
            family.validateResponse(response[row], weights[row]);
        }

        GlmResult initial = Glm.fit(response, fixedEffects, rows, columns,
            family, weights, offsets, options.initialGlmOptions(), backendPolicy);
        double[] linearPredictor = initial.linearPredictor();
        double[] previousVariances = options.remlOptions().initialVariances();
        RemlResult workingFit = null;
        double[] randomPredictor = new double[rows];
        double[] fittedMeans = initial.fittedMeans();
        double conditionalDeviance = deviance(
            response, fittedMeans, weights, family);
        boolean converged = false;
        String message = "maximum PQL iterations reached";
        int iterations = 0;

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            for (int iteration = 1;
                    iteration <= options.maximumIterations(); iteration++) {
                iterations = iteration;
                WorkingModel working = workingModel(response, linearPredictor,
                    weights, offsets, family);
                RemlOptions innerOptions = withInitialVariances(
                    options.remlOptions(), previousVariances);
                workingFit = Reml.fitWithKnownCovariance(
                    working.response(), fixedEffects, rows, columns,
                    randomComponents, working.residualCovariance(),
                    innerOptions, context);

                double[] variances = workingFit.varianceComponents();
                double[] fixedPredictor = MatrixOps.multiply(backend,
                    fixedEffects, rows, columns, workingFit.fixedEffects());
                double[] randomCovariance = covariance(
                    randomComponents, variances, rows);
                double[] totalCovariance = randomCovariance.clone();
                for (int index = 0; index < totalCovariance.length; index++) {
                    totalCovariance[index] += working.residualCovariance()[index];
                }
                CholeskyFactor factor = backend.dpotrf(totalCovariance, rows);
                double[] workingResidual = new double[rows];
                for (int row = 0; row < rows; row++) {
                    workingResidual[row] = working.response()[row]
                        - fixedPredictor[row];
                }
                double[] projected = factor.solve(workingResidual);
                double[] candidateRandom = MatrixOps.multiply(
                    backend, randomCovariance, rows, rows, projected);
                double[] candidatePredictor = new double[rows];
                double[] candidateMeans = new double[rows];
                for (int row = 0; row < rows; row++) {
                    candidatePredictor[row] = offsets[row] + fixedPredictor[row]
                        + candidateRandom[row];
                    candidateMeans[row] = family.inverseLink(candidatePredictor[row]);
                }

                double predictorChange = relativeMaximumChange(
                    linearPredictor, candidatePredictor);
                double varianceChange = previousVariances == null
                    ? Double.POSITIVE_INFINITY
                    : relativeMaximumChange(previousVariances, variances);
                linearPredictor = candidatePredictor;
                randomPredictor = candidateRandom;
                fittedMeans = candidateMeans;
                previousVariances = variances;
                conditionalDeviance = deviance(
                    response, fittedMeans, weights, family);

                if (predictorChange <= options.relativeTolerance()
                        && varianceChange <= options.relativeTolerance()
                        && workingFit.converged()) {
                    converged = true;
                    message = "predictor and variance tolerances reached";
                    break;
                }
            }

            if (workingFit == null) {
                throw new IllegalStateException("PQL did not execute a working-model fit");
            }
            if (!workingFit.converged() && !converged) {
                message = "final REML working model did not converge: "
                    + workingFit.convergenceMessage();
            }

            List<String> names = new ArrayList<>(randomComponents.size());
            for (VarianceComponent component : randomComponents) {
                names.add(component.name());
            }
            return new GlmmPqlResult(
                family.name(), names, workingFit.varianceComponents(),
                workingFit.fixedEffects(), workingFit.fixedEffectCovariance(),
                workingFit.fixedEffectInferenceCovariance(),
                workingFit.fixedEffectStandardErrors(),
                workingFit.associationStatistics(), randomPredictor,
                linearPredictor, fittedMeans, conditionalDeviance,
                workingFit.restrictedLogLikelihood(), iterations,
                converged, message, workingFit.backend());
        }
    }

    private static WorkingModel workingModel(
            double[] response, double[] predictor,
            double[] priorWeights, double[] offset, GlmFamily family) {
        int rows = response.length;
        double[] workingResponse = new double[rows];
        double[] covariance = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            double mean = family.inverseLink(predictor[row]);
            double derivative = family.meanDerivative(predictor[row]);
            double variance = family.variance(mean);
            if (!Double.isFinite(derivative) || derivative == 0.0
                    || !Double.isFinite(variance) || variance <= 0.0) {
                throw new IllegalArgumentException(
                    "family produced invalid PQL derivative or variance");
            }
            double weight = clamp(priorWeights[row]
                * derivative * derivative / variance,
                MINIMUM_WORKING_WEIGHT, MAXIMUM_WORKING_WEIGHT);
            workingResponse[row] = predictor[row]
                + (response[row] - mean) / derivative - offset[row];
            covariance[row * rows + row] = 1.0 / weight;
        }
        return new WorkingModel(workingResponse, covariance);
    }

    private static RemlOptions withInitialVariances(
            RemlOptions source, double[] variances) {
        RemlOptions.Builder builder = RemlOptions.builder()
            .maximumIterations(source.maximumIterations())
            .relativeTolerance(source.relativeTolerance())
            .scoreTolerance(source.scoreTolerance())
            .varianceBounds(source.minimumVariance(), source.maximumVariance())
            .maximumLogVarianceStep(source.maximumLogVarianceStep())
            .degreesOfFreedomMethod(source.degreesOfFreedomMethod())
            .varianceEstimation(source.varianceEstimation());
        if (variances != null) {
            builder.initialVariances(variances);
        }
        return builder.build();
    }

    private static double[] covariance(
            List<VarianceComponent> components,
            double[] variances,
            int dimension) {
        double[] result = new double[dimension * dimension];
        for (int component = 0; component < components.size(); component++) {
            double[] basis = components.get(component).covariance();
            for (int index = 0; index < result.length; index++) {
                result[index] += variances[component] * basis[index];
            }
        }
        return result;
    }

    private static double deviance(
            double[] response, double[] means,
            double[] weights, GlmFamily family) {
        double result = 0.0;
        for (int row = 0; row < response.length; row++) {
            result += weights[row]
                * family.unitDeviance(response[row], means[row]);
        }
        return result;
    }

    private static double relativeMaximumChange(
            double[] previous, double[] current) {
        double maximum = 0.0;
        for (int index = 0; index < previous.length; index++) {
            maximum = Math.max(maximum,
                Math.abs(current[index] - previous[index])
                    / (1.0 + Math.abs(previous[index])));
        }
        return maximum;
    }

    private static void validateComponents(
            List<VarianceComponent> components, int dimension) {
        if (components == null || components.isEmpty()) {
            throw new IllegalArgumentException(
                "at least one random covariance component is required");
        }
        for (VarianceComponent component : components) {
            if (component == null || component.dimension() != dimension) {
                throw new IllegalArgumentException(
                    "random components must match the response dimension");
            }
        }
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

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    private record WorkingModel(
            double[] response,
            double[] residualCovariance) {
    }
}
