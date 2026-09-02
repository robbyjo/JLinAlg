/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glmm.GlmmPqlResult;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.VarianceComponent;

/** Reconstructs the final PQL working model for modes, EDF, and prediction. */
final class PqlWorkingState {
    private static final double MINIMUM_WEIGHT = 1e-12;
    private static final double MAXIMUM_WEIGHT = 1e150;

    private PqlWorkingState() { }

    static State reconstruct(
            double[] response,
            double[] fixed,
            int rows,
            int columns,
            GlmFamily family,
            double[] priorWeights,
            double[] offset,
            List<VarianceComponent> components,
            GlmmPqlResult fitted,
            ComputeBackend backend) {
        double[] predictor = fitted.linearPredictor();
        double[] means = fitted.fittedMeans();
        double[] workingResponse = new double[rows];
        double[] covariance = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            double derivative = family.meanDerivative(predictor[row]);
            double variance = family.variance(means[row]);
            double weight = clamp(priorWeights[row]
                * derivative * derivative / variance,
                MINIMUM_WEIGHT, MAXIMUM_WEIGHT);
            workingResponse[row] = predictor[row]
                + (response[row] - means[row]) / derivative - offset[row];
            covariance[row * rows + row] = 1.0 / weight;
        }
        double[] variances = fitted.varianceComponents();
        for (int component = 0; component < components.size(); component++) {
            double[] basis = components.get(component).covariance();
            for (int index = 0; index < covariance.length; index++) {
                covariance[index] += variances[component] * basis[index];
            }
        }
        CholeskyFactor factor = backend.dpotrf(covariance, rows);
        double[] fixedPredictor = MatrixOps.multiply(
            backend, fixed, rows, columns, fitted.beta());
        double[] projected = factor.solve(
            MatrixOps.subtract(workingResponse, fixedPredictor));
        double[] inverse = factor.solve(MatrixOps.identity(rows), rows);
        double[] inverseFixed = factor.solve(fixed, columns);
        double[] temporary = MatrixOps.multiply(backend,
            inverseFixed, rows, columns,
            fitted.fixedEffectCovariance(), columns);
        double[] correction = new double[rows * rows];
        backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
            rows, rows, columns, 1.0,
            temporary, inverseFixed, 0.0, correction);
        return new State(projected, MatrixOps.subtract(inverse, correction));
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    record State(double[] projectedResidual, double[] projection) { }
}
