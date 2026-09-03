/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPql;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/**
 * Retained first-order PQL efficient-score projection for non-Gaussian set
 * tests. Null fitting and the dense projection are performed once.
 */
public final class PqlSetTestNullModel implements SetTestScoreNullModel {
    private static final double MINIMUM_WORKING_WEIGHT = 1e-12;
    private static final double MAXIMUM_WORKING_WEIGHT = 1e150;

    private final double[] projection;
    private final double[] projectedResponse;
    private final int observations;
    private final GlmmPqlResult nullModel;
    private final BackendPolicy backendPolicy;

    private PqlSetTestNullModel(
            double[] projection, double[] projectedResponse,
            int observations, GlmmPqlResult nullModel,
            BackendPolicy backendPolicy) {
        this.projection = projection;
        this.projectedResponse = projectedResponse;
        this.observations = observations;
        this.nullModel = nullModel;
        this.backendPolicy = backendPolicy;
    }

    public static PqlSetTestNullModel prepare(
            double[] response, double[][] fixedEffects, GlmFamily family,
            List<VarianceComponent> randomComponents,
            GlmmPqlOptions options, BackendPolicy backendPolicy) {
        if (response == null)
            throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(fixedEffects, response.length);
        return prepare(response, fixed, response.length,
            fixedEffects[0].length, family, randomComponents, options,
            backendPolicy);
    }

    public static PqlSetTestNullModel prepare(
            double[] response, double[] fixedEffects, int rows, int columns,
            GlmFamily family, List<VarianceComponent> randomComponents,
            GlmmPqlOptions options, BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        if (family == null || randomComponents == null
                || randomComponents.isEmpty() || options == null
                || backendPolicy == null)
            throw new IllegalArgumentException(
                "family, random components, and options are required");
        GlmmPqlResult fit = GlmmPql.fit(response, fixedEffects, rows, columns,
            family, randomComponents, null, null, options, backendPolicy);
        if (!fit.converged())
            throw new IllegalArgumentException(
                "PQL set-test null did not converge: "
                    + fit.convergenceMessage());

        double[] predictor = fit.linearPredictor();
        double[] means = fit.fittedMeans();
        double[] workingResponse = new double[rows];
        double[] covariance = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            double derivative = family.meanDerivative(predictor[row]);
            double variance = family.variance(means[row]);
            if (!Double.isFinite(derivative) || derivative == 0.0
                    || !(variance > 0.0) || !Double.isFinite(variance))
                throw new IllegalArgumentException(
                    "family produced invalid final PQL working weights");
            double weight = clamp(derivative * derivative / variance,
                MINIMUM_WORKING_WEIGHT, MAXIMUM_WORKING_WEIGHT);
            workingResponse[row] = predictor[row]
                + (response[row] - means[row]) / derivative;
            covariance[row * rows + row] = 1.0 / weight;
        }
        double[] variances = fit.varianceComponents();
        for (int component = 0; component < randomComponents.size(); component++) {
            double[] basis = randomComponents.get(component).covariance();
            for (int index = 0; index < covariance.length; index++)
                covariance[index] += variances[component] * basis[index];
        }

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            CholeskyFactor factor = backend.dpotrf(covariance, rows);
            double[] inverse = factor.solve(MatrixOps.identity(rows), rows);
            double[] inverseFixed = factor.solve(fixedEffects, columns);
            double[] temporary = MatrixOps.multiply(backend,
                inverseFixed, rows, columns,
                fit.fixedEffectCovariance(), columns);
            double[] correction = new double[rows * rows];
            backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                rows, rows, columns, 1.0,
                temporary, inverseFixed, 0.0, correction);
            double[] projection = MatrixOps.subtract(inverse, correction);
            double[] fixedPredictor = MatrixOps.multiply(backend,
                fixedEffects, rows, columns, fit.fixedEffects());
            for (int row = 0; row < rows; row++)
                workingResponse[row] -= fixedPredictor[row];
            double[] projectedResponse = factor.solve(workingResponse);
            return new PqlSetTestNullModel(projection, projectedResponse,
                rows, fit, backendPolicy);
        }
    }

    @Override public int observations() { return observations; }

    @Override
    public SetTestScoreState score(double[][] variantRows) {
        if (variantRows == null || variantRows.length == 0)
            throw new IllegalArgumentException("variant rows are required");
        int variants = variantRows.length;
        double[] sampleByVariant = new double[observations * variants];
        for (int variant = 0; variant < variants; variant++) {
            if (variantRows[variant] == null
                    || variantRows[variant].length != observations)
                throw new IllegalArgumentException(
                    "variant rows must match null-model observations");
            MatrixOps.requireFinite(variantRows[variant], "variant row");
            for (int sample = 0; sample < observations; sample++)
                sampleByVariant[sample * variants + variant] =
                    variantRows[variant][sample];
        }
        double[] projected;
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            projected = MatrixOps.multiply(context.backend(), projection,
                observations, observations, sampleByVariant, variants);
        }
        double[] scores = new double[variants];
        double[] information = new double[variants * variants];
        for (int left = 0; left < variants; left++) {
            for (int sample = 0; sample < observations; sample++)
                scores[left] += variantRows[left][sample]
                    * projectedResponse[sample];
            for (int right = 0; right <= left; right++) {
                double value = 0;
                for (int sample = 0; sample < observations; sample++)
                    value += variantRows[left][sample]
                        * projected[sample * variants + right];
                information[left * variants + right] = value;
                information[right * variants + left] = value;
            }
        }
        return new SetTestScoreState(scores, information, variants);
    }

    public GlmmPqlResult nullModel() { return nullModel; }
    @Override public BackendPolicy backendPolicy() { return backendPolicy; }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }
}
