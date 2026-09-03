/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixDiagonal;
import jdistlib.accelerator.MatrixSide;
import jdistlib.accelerator.MatrixTranspose;
import jdistlib.accelerator.MatrixTriangle;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Retained efficient-score projection for a caller-supplied null covariance. */
public final class KnownCovarianceSetTestNullModel
        implements SetTestScoreNullModel, AutoCloseable {
    private final BackendContext context;
    private final double[] projection;
    private final double[] projectedResponse;
    private final int observations;
    private final BackendPolicy backendPolicy;

    private KnownCovarianceSetTestNullModel(
            BackendContext context, double[] projection,
            double[] projectedResponse, int observations,
            BackendPolicy backendPolicy) {
        this.context = context;
        this.projection = projection;
        this.projectedResponse = projectedResponse;
        this.observations = observations;
        this.backendPolicy = backendPolicy;
    }

    /**
     * Prepares P and Py. For Gaussian models, {@code workingResponse} is y and
     * {@code covariance} is the fitted marginal covariance. For a GLMM it is
     * the final pseudo-response and PQL working covariance.
     */
    public static KnownCovarianceSetTestNullModel prepare(
            double[] workingResponse, double[][] fixedEffects,
            double[] covariance, BackendPolicy backendPolicy) {
        if (workingResponse == null)
            throw new IllegalArgumentException("working response is required");
        double[] fixed = MatrixOps.rowMajor(
            fixedEffects, workingResponse.length);
        return prepare(workingResponse, fixed, workingResponse.length,
            fixedEffects[0].length, covariance, backendPolicy);
    }

    public static KnownCovarianceSetTestNullModel prepare(
            double[] workingResponse, double[] fixedEffects,
            int rows, int columns, double[] covariance,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(
            workingResponse, fixedEffects, rows, columns);
        if (covariance == null || covariance.length != rows * rows
                || backendPolicy == null)
            throw new IllegalArgumentException(
                "covariance and backend policy are required");
        MatrixOps.requireFinite(covariance, "covariance");
        BackendContext context = BackendContext.select(backendPolicy);
        try {
            ComputeBackend backend = context.backend();
            CholeskyFactor factor = backend.dpotrf(covariance, rows);
            double[] lower = factor.lower();
            double[] inverse = MatrixOps.identity(rows);
            solve(backend, lower, rows, inverse, rows);
            double[] inverseFixed = fixedEffects.clone();
            solve(backend, lower, rows, inverseFixed, columns);
            double[] information = MatrixOps.transposeMultiply(backend,
                fixedEffects, rows, columns, inverseFixed, columns);
            CholeskyFactor informationFactor = backend.dpotrf(
                information, columns);
            double[] fixedCovariance = informationFactor.solve(
                MatrixOps.identity(columns), columns);
            double[] temporary = MatrixOps.multiply(backend,
                inverseFixed, rows, columns, fixedCovariance, columns);
            double[] correction = new double[rows * rows];
            backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
                rows, rows, columns, 1.0,
                temporary, inverseFixed, 0.0, correction);
            double[] projection = MatrixOps.subtract(inverse, correction);
            double[] projectedResponse = MatrixOps.multiply(backend,
                projection, rows, rows, workingResponse);
            return new KnownCovarianceSetTestNullModel(
                context, projection, projectedResponse, rows, backendPolicy);
        } catch (RuntimeException | Error failure) {
            context.close();
            throw failure;
        }
    }

    @Override public int observations() { return observations; }

    @Override
    public synchronized SetTestScoreState score(double[][] variantRows) {
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
        ComputeBackend backend = context.backend();
        double[] projected = MatrixOps.multiply(backend, projection,
            observations, observations, sampleByVariant, variants);
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

    public BackendPolicy backendPolicy() { return backendPolicy; }

    @Override public ComputeBackend computeBackend() {
        return context.backend();
    }

    @Override public void close() {
        context.close();
    }

    private static void solve(
            ComputeBackend backend, double[] lower, int dimension,
            double[] rightHandSide, int columns) {
        backend.dtrsm(MatrixSide.LEFT, MatrixTriangle.LOWER,
            MatrixTranspose.NONE, MatrixDiagonal.NON_UNIT,
            dimension, columns, 1.0, lower, rightHandSide);
        backend.dtrsm(MatrixSide.LEFT, MatrixTriangle.LOWER,
            MatrixTranspose.TRANSPOSE, MatrixDiagonal.NON_UNIT,
            dimension, columns, 1.0, lower, rightHandSide);
    }
}
