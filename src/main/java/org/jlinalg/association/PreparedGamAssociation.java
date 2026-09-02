/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.gam.GaussianSmoothSelectionResult;
import org.jlinalg.gam.PenalizedPredictor;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;

/** Retained penalized-null projection for fast GAM-adjusted association scans. */
public final class PreparedGamAssociation implements AutoCloseable {
    private final int observations;
    private final int columns;
    private final double[] design;
    private final double[] residuals;
    private final double residualVariance;
    private final double degreesOfFreedom;
    private final AssociationEngineOptions options;
    private final BackendContext context;
    private final CholeskyFactor penalizedFactor;
    private boolean closed;

    /** Retains the fitted penalty matrix and factorizes its normal equations once. */
    public PreparedGamAssociation(
            GaussianSmoothSelectionResult nullModel,
            AssociationEngineOptions options) {
        if (nullModel == null || options == null) {
            throw new IllegalArgumentException("null model and engine options are required");
        }
        this.options = options;
        PenalizedPredictor predictor = nullModel.predictor();
        this.observations = predictor.observations();
        this.columns = predictor.columns();
        this.design = predictor.design();
        this.residuals = nullModel.residuals();
        this.residualVariance = nullModel.residualVariance();
        this.degreesOfFreedom = observations
            - nullModel.effectiveDegreesOfFreedom() - 1.0;
        if (!(degreesOfFreedom > 0.0)) {
            throw new IllegalArgumentException(
                "GAM association scan requires positive residual degrees of freedom");
        }
        context = BackendContext.select(options.backendPolicy());
        ComputeBackend backend = context.backend();
        double[] cross = MatrixOps.transposeMultiply(
            backend, design, observations, columns, design, columns);
        double[] penalty = predictor.penaltyDiagonal();
        for (int column = 0; column < columns; column++) {
            cross[column * columns + column] += penalty[column];
        }
        penalizedFactor = backend.dpotrf(cross, columns);
    }

    /** Scans columns of an observation-by-variable predictor matrix. */
    public AssociationBatchResult scan(
            double[][] predictors, List<String> names) {
        requireOpen();
        if (predictors == null || predictors.length != observations
                || predictors[0] == null || predictors[0].length == 0) {
            throw new IllegalArgumentException(
                "predictors must have one nonempty row per observation");
        }
        int predictorCount = predictors[0].length;
        if (names == null || names.size() != predictorCount) {
            throw new IllegalArgumentException("one name is required per predictor");
        }
        double[] changing = MatrixOps.rowMajorUnchecked(predictors, observations);
        imputeOrValidate(changing, predictorCount);
        long started = System.nanoTime();
        ComputeBackend backend = context.backend();
        double[] cross = MatrixOps.transposeMultiply(backend,
            design, observations, columns, changing, predictorCount);
        double[] solved = penalizedFactor.solve(cross, predictorCount);
        double[] beta = new double[predictorCount];
        double[] standardErrors = new double[predictorCount];
        java.util.Arrays.fill(beta, Double.NaN);
        java.util.Arrays.fill(standardErrors, Double.NaN);
        List<AssociationFailure> failures =
            java.util.Collections.synchronizedList(new ArrayList<>());
        Runnable work = () -> IntStream.range(0, predictorCount).parallel()
            .forEach(index -> calculate(index, names.get(index), changing,
                predictorCount, cross, solved, beta, standardErrors, failures));
        if (options.parallelism() == 1 || predictorCount == 1) {
            IntStream.range(0, predictorCount).forEach(index -> calculate(index,
                names.get(index), changing, predictorCount, cross, solved,
                beta, standardErrors, failures));
        } else {
            ForkJoinPool pool = new ForkJoinPool(options.parallelism());
            try {
                pool.submit(work).join();
            } finally {
                pool.shutdown();
            }
        }
        AssociationStatistics statistics = AssociationStatistics.studentT(
            beta, standardErrors, degreesOfFreedom,
            DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
        return new AssociationBatchResult(List.copyOf(names),
            statistics.beta(), statistics.standardErrors(),
            statistics.statistics(), statistics.pValues(),
            statistics.degreesOfFreedom(),
            statistics.negativeLog10PValues(), List.copyOf(failures),
            -1, options.parallelism(), System.nanoTime() - started);
    }

    private void calculate(
            int index,
            String name,
            double[] predictors,
            int predictorCount,
            double[] cross,
            double[] solved,
            double[] beta,
            double[] standardErrors,
            List<AssociationFailure> failures) {
        double numerator = 0.0;
        double rawInformation = 0.0;
        for (int row = 0; row < observations; row++) {
            double value = predictors[row * predictorCount + index];
            numerator += value * residuals[row];
            rawInformation += value * value;
        }
        double removed = 0.0;
        for (int column = 0; column < columns; column++) {
            removed += cross[column * predictorCount + index]
                * solved[column * predictorCount + index];
        }
        double information = rawInformation - removed;
        if (!(information > 1e-12 * Math.max(1.0, rawInformation))
                || !Double.isFinite(information)) {
            String message = "predictor has no information after GAM adjustment";
            if (options.failurePolicy() == AssociationFailurePolicy.FAIL_FAST) {
                throw new IllegalArgumentException(name + ": " + message);
            }
            failures.add(new AssociationFailure(index, name,
                IllegalArgumentException.class.getName(), message));
            return;
        }
        beta[index] = numerator / information;
        standardErrors[index] = Math.sqrt(residualVariance / information);
    }

    private void imputeOrValidate(double[] values, int predictorCount) {
        for (int column = 0; column < predictorCount; column++) {
            double sum = 0.0;
            int count = 0;
            for (int row = 0; row < observations; row++) {
                double value = values[row * predictorCount + column];
                if (Double.isFinite(value)) {
                    sum += value;
                    count++;
                } else if (options.predictorMissingPolicy()
                        == VariableMissingPolicy.ERROR) {
                    throw new IllegalArgumentException(
                        "predictor " + column + " contains a missing value");
                }
            }
            if (count == 0) {
                throw new IllegalArgumentException(
                    "predictor " + column + " is entirely missing");
            }
            double mean = sum / count;
            for (int row = 0; row < observations; row++) {
                int index = row * predictorCount + column;
                if (!Double.isFinite(values[index])) values[index] = mean;
            }
        }
    }

    public int observations() { return observations; }
    public double degreesOfFreedom() { return degreesOfFreedom; }

    private void requireOpen() {
        if (closed) throw new IllegalStateException("prepared GAM scan is closed");
    }

    @Override public void close() {
        if (!closed) {
            closed = true;
            context.close();
        }
    }
}
