/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.ols.RankDeficiencyStrategy;

/** Covariate-residualized, block-matrix OLS association scans. */
public final class FastOlsAssociation {
    private FastOlsAssociation() { }

    public static AssociationBatchResult scanPredictors(
            double[] response,
            double[][] covariates,
            double[][] predictors,
            List<String> predictorNames,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) {
        if (response == null) throw new IllegalArgumentException("response is required");
        double[] fixed = MatrixOps.rowMajor(covariates, response.length);
        double[] changing = MatrixOps.rowMajorUnchecked(predictors, response.length);
        return scanPredictors(response, fixed, response.length, covariates[0].length,
            changing, predictors[0].length, predictorNames,
            weights, offset, olsOptions, engineOptions);
    }

    /** Fast weighted predictor scan from contiguous row-major matrices. */
    public static AssociationBatchResult scanPredictors(
            double[] response,
            double[] covariates,
            int rows,
            int covariateCount,
            double[] predictors,
            int predictorCount,
            List<String> predictorNames,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) {
        MatrixOps.validateModelData(response, covariates, rows, covariateCount);
        validate(predictors, rows, predictorCount, olsOptions, engineOptions);
        List<String> names = names(predictorNames, predictorCount, "variable");
        long started = System.nanoTime();
        WeightedBase base = prepareBase(response, covariates, rows,
            covariateCount, weights, offset, olsOptions, engineOptions);
        int degreesOfFreedom = rows - base.rank() - 1;
        if (degreesOfFreedom < 1)
            throw new IllegalArgumentException("predictor scan requires residual degrees of freedom");
        double[] beta = nan(predictorCount);
        double[] standardErrors = nan(predictorCount);
        ConcurrentLinkedQueue<AssociationFailure> failures = new ConcurrentLinkedQueue<>();
        int chunks = (predictorCount + engineOptions.chunkSize() - 1)
            / engineOptions.chunkSize();
        execute(chunks, engineOptions.parallelism(), chunk -> {
            int first = chunk * engineOptions.chunkSize();
            int count = Math.min(engineOptions.chunkSize(), predictorCount - first);
            try (BackendContext context = BackendContext.select(
                    engineOptions.backendPolicy())) {
                ComputeBackend backend = context.backend();
                PreparedPredictors prepared = predictorBlock(
                    predictors, rows, predictorCount,
                    first, count, base.squareRootWeights(),
                    engineOptions.predictorMissingPolicy(), names, failures,
                    engineOptions.failurePolicy());
                double[] block = prepared.values();
                double[] fixedCross = MatrixOps.transposeMultiply(backend,
                    base.weightedDesign(), rows, covariateCount, block, count);
                double[] projectionCoefficients = MatrixOps.multiply(backend,
                    base.informationInverse(), covariateCount, covariateCount,
                    fixedCross, count);
                double[] fixedProjection = MatrixOps.multiply(backend,
                    base.weightedDesign(), rows, covariateCount,
                    projectionCoefficients, count);
                for (int index = 0; index < block.length; index++)
                    block[index] -= fixedProjection[index];
                for (int marker = 0; marker < count; marker++) {
                    int destination = first + marker;
                    if (!prepared.valid()[marker]) continue;
                    double numerator = 0.0;
                    double information = 0.0;
                    for (int row = 0; row < rows; row++) {
                        double value = block[row * count + marker];
                        numerator += value * base.responseResidual()[row];
                        information += value * value;
                    }
                    if (!(information > 1e-14) || !Double.isFinite(information)) {
                        failures.add(new AssociationFailure(destination,
                            names.get(destination), "NonEstimablePredictor",
                            "predictor is constant or collinear with fixed covariates"));
                        if (engineOptions.failurePolicy() == AssociationFailurePolicy.FAIL_FAST)
                            throw new IllegalArgumentException(
                                "predictor is not estimable: " + names.get(destination));
                        continue;
                    }
                    beta[destination] = numerator / information;
                    double rss = Math.max(0.0, base.residualSumSquares()
                        - numerator * numerator / information);
                    standardErrors[destination] = Math.sqrt(
                        rss / degreesOfFreedom / information);
                }
            }
        });
        AssociationStatistics statistics = AssociationStatistics.studentT(
            beta, standardErrors, degreesOfFreedom, DegreesOfFreedomMethod.RESIDUAL);
        return result(names, statistics, failures, covariateCount,
            Math.min(engineOptions.parallelism(), chunks),
            System.nanoTime() - started);
    }

    public static AssociationBatchResult scanResponses(
            double[][] responses,
            double[][] design,
            int coefficientIndex,
            List<String> responseNames,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) {
        if (responses == null || responses.length == 0)
            throw new IllegalArgumentException("responses are required");
        double[] values = MatrixOps.rowMajorUnchecked(responses, responses.length);
        double[] fixed = MatrixOps.rowMajor(design, responses.length);
        return scanResponses(values, responses.length, responses[0].length,
            fixed, design[0].length, coefficientIndex, responseNames,
            weights, offset, olsOptions, engineOptions);
    }

    /** Reuses one weighted fixed-design information inverse across response traits. */
    public static AssociationBatchResult scanResponses(
            double[] responses,
            int rows,
            int responseCount,
            double[] design,
            int columns,
            int coefficientIndex,
            List<String> responseNames,
            double[] weights,
            double[] offset,
            OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) {
        if (responses == null || responses.length != rows * responseCount
                || responseCount < 1 || coefficientIndex < 0
                || coefficientIndex >= columns || olsOptions == null
                || engineOptions == null) {
            throw new IllegalArgumentException("response scan dimensions or options are invalid");
        }
        MatrixOps.validateModelData(new double[rows], design, rows, columns);
        double[] squareRootWeights = weights(weights, rows);
        double[] offsets = offsets(offset, rows);
        double[] weightedDesign = weightDesign(design, rows, columns, squareRootWeights);
        LeastSquaresSolver.Solution factor;
        try (BackendContext context = BackendContext.select(engineOptions.backendPolicy())) {
            factor = LeastSquaresSolver.solve(weightedDesign, new double[rows],
                rows, columns,
                olsOptions.rankDeficiencyStrategy() == RankDeficiencyStrategy.MINIMUM_NORM,
                context.backend());
        }
        int degreesOfFreedom = rows - factor.rank();
        if (degreesOfFreedom < 1)
            throw new IllegalArgumentException("response scan requires residual degrees of freedom");
        double inverseDiagonal = factor.unscaledCovariance()[
            coefficientIndex * columns + coefficientIndex];
        List<String> names = names(responseNames, responseCount, "response");
        double[] beta = nan(responseCount);
        double[] standardErrors = nan(responseCount);
        ConcurrentLinkedQueue<AssociationFailure> failures = new ConcurrentLinkedQueue<>();
        long started = System.nanoTime();
        int chunks = (responseCount + engineOptions.chunkSize() - 1)
            / engineOptions.chunkSize();
        execute(chunks, engineOptions.parallelism(), chunk -> {
            double[] weightedResponse = new double[rows];
            try (BackendContext context = BackendContext.select(engineOptions.backendPolicy())) {
                ComputeBackend backend = context.backend();
                for (int trait = chunk * engineOptions.chunkSize();
                        trait < Math.min(responseCount,
                            (chunk + 1) * engineOptions.chunkSize()); trait++) {
                    try {
                        for (int row = 0; row < rows; row++) {
                            double value = responses[row * responseCount + trait];
                            if (!Double.isFinite(value))
                                throw new IllegalArgumentException("response contains a non-finite value");
                            weightedResponse[row] = squareRootWeights[row]
                                * (value - offsets[row]);
                        }
                        double[] cross = new double[columns];
                        backend.dgemv(jdistlib.accelerator.MatrixTranspose.TRANSPOSE,
                            rows, columns, 1.0, weightedDesign, weightedResponse,
                            0.0, cross);
                        double[] coefficients = MatrixOps.multiply(backend,
                            factor.unscaledCovariance(), columns, columns, cross);
                        double[] fitted = MatrixOps.multiply(backend,
                            weightedDesign, rows, columns, coefficients);
                        double rss = 0.0;
                        for (int row = 0; row < rows; row++) {
                            double residual = weightedResponse[row] - fitted[row];
                            rss += residual * residual;
                        }
                        beta[trait] = coefficients[coefficientIndex];
                        standardErrors[trait] = Math.sqrt(
                            Math.max(0.0, rss / degreesOfFreedom * inverseDiagonal));
                    } catch (RuntimeException exception) {
                        if (engineOptions.failurePolicy() == AssociationFailurePolicy.FAIL_FAST)
                            throw exception;
                        failures.add(failure(trait, names.get(trait), exception));
                    }
                }
            }
        });
        AssociationStatistics statistics = AssociationStatistics.studentT(
            beta, standardErrors, degreesOfFreedom, DegreesOfFreedomMethod.RESIDUAL);
        return result(names, statistics, failures, coefficientIndex,
            Math.min(engineOptions.parallelism(), chunks),
            System.nanoTime() - started);
    }

    private static WeightedBase prepareBase(
            double[] response, double[] design, int rows, int columns,
            double[] weights, double[] offset, OlsOptions olsOptions,
            AssociationEngineOptions engineOptions) {
        double[] squareRootWeights = weights(weights, rows);
        double[] offsets = offsets(offset, rows);
        double[] weightedDesign = weightDesign(design, rows, columns, squareRootWeights);
        double[] weightedResponse = new double[rows];
        for (int row = 0; row < rows; row++)
            weightedResponse[row] = squareRootWeights[row] * (response[row] - offsets[row]);
        try (BackendContext context = BackendContext.select(engineOptions.backendPolicy())) {
            ComputeBackend backend = context.backend();
            LeastSquaresSolver.Solution solution = LeastSquaresSolver.solve(
                weightedDesign, weightedResponse, rows, columns,
                olsOptions.rankDeficiencyStrategy() == RankDeficiencyStrategy.MINIMUM_NORM,
                backend);
            double[] fitted = MatrixOps.multiply(
                backend, weightedDesign, rows, columns, solution.coefficients());
            double[] residual = MatrixOps.subtract(weightedResponse, fitted);
            double rss = backend.ddot(rows, residual, 0, 1, residual, 0, 1);
            return new WeightedBase(weightedDesign,
                solution.unscaledCovariance(), residual, squareRootWeights,
                Math.max(0.0, rss), solution.rank());
        }
    }

    private static PreparedPredictors predictorBlock(
            double[] predictors, int rows, int total, int first, int count,
            double[] squareRootWeights, VariableMissingPolicy policy,
            List<String> names,
            ConcurrentLinkedQueue<AssociationFailure> failures,
            AssociationFailurePolicy failurePolicy) {
        double[] result = new double[rows * count];
        boolean[] valid = new boolean[count];
        Arrays.fill(valid, true);
        for (int marker = 0; marker < count; marker++) {
            int source = first + marker;
            double sum = 0.0;
            int finite = 0;
            int firstMissing = -1;
            for (int row = 0; row < rows; row++) {
                double value = predictors[row * total + source];
                if (Double.isFinite(value)) { sum += value; finite++; }
                else if (firstMissing < 0) firstMissing = row;
            }
            String message = null;
            if (finite == 0) message = "predictor has no finite values";
            else if (firstMissing >= 0 && policy == VariableMissingPolicy.ERROR)
                message = "non-finite predictor at row " + firstMissing;
            if (message != null) {
                if (failurePolicy == AssociationFailurePolicy.FAIL_FAST)
                    throw new IllegalArgumentException(
                        names.get(source) + ": " + message);
                valid[marker] = false;
                failures.add(new AssociationFailure(source, names.get(source),
                    "IllegalArgumentException", message));
                continue;
            }
            double mean = sum / finite;
            for (int row = 0; row < rows; row++) {
                double value = predictors[row * total + source];
                result[row * count + marker] = squareRootWeights[row]
                    * (Double.isFinite(value) ? value : mean);
            }
        }
        return new PreparedPredictors(result, valid);
    }

    private static double[] weightDesign(
            double[] design, int rows, int columns, double[] squareRootWeights) {
        double[] result = new double[design.length];
        for (int row = 0; row < rows; row++)
            for (int column = 0; column < columns; column++)
                result[row * columns + column] = squareRootWeights[row]
                    * design[row * columns + column];
        return result;
    }

    private static double[] weights(double[] weights, int rows) {
        double[] result = new double[rows];
        if (weights == null) Arrays.fill(result, 1.0);
        else {
            if (weights.length != rows) throw new IllegalArgumentException("weight length must equal rows");
            for (int row = 0; row < rows; row++) {
                if (!(weights[row] > 0.0) || !Double.isFinite(weights[row]))
                    throw new IllegalArgumentException("weights must be finite and positive");
                result[row] = Math.sqrt(weights[row]);
            }
        }
        return result;
    }

    private static double[] offsets(double[] offset, int rows) {
        if (offset == null) return new double[rows];
        if (offset.length != rows) throw new IllegalArgumentException("offset length must equal rows");
        return MatrixOps.finiteCopy(offset, "offset");
    }

    private static void validate(
            double[] predictors, int rows, int count,
            OlsOptions ols, AssociationEngineOptions engine) {
        if (predictors == null || count < 1 || predictors.length != rows * count
                || ols == null || engine == null)
            throw new IllegalArgumentException("predictors or options are invalid");
    }

    private static void execute(int chunks, int parallelism, Chunk operation) {
        ForkJoinPool pool = new ForkJoinPool(Math.min(chunks, parallelism));
        try {
            pool.submit(() -> IntStream.range(0, chunks).parallel()
                .forEach(operation::run)).get();
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("association scan was interrupted", exception);
        } catch (ExecutionException exception) {
            Throwable cause = exception.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            throw new IllegalStateException("association scan failed", cause);
        } finally { pool.shutdownNow(); }
    }

    private static AssociationBatchResult result(
            List<String> names, AssociationStatistics statistics,
            ConcurrentLinkedQueue<AssociationFailure> failures,
            int coefficient, int parallelism, long elapsed) {
        List<AssociationFailure> ordered = failures.stream()
            .sorted(Comparator.comparingInt(AssociationFailure::index)).toList();
        return new AssociationBatchResult(names, statistics.beta(),
            statistics.standardErrors(), statistics.statistics(),
            statistics.pValues(), statistics.degreesOfFreedom(),
            statistics.negativeLog10PValues(), ordered,
            coefficient, parallelism, elapsed);
    }

    private static AssociationFailure failure(
            int index, String name, RuntimeException exception) {
        return new AssociationFailure(index, name,
            exception.getClass().getSimpleName(),
            exception.getMessage() == null ? "" : exception.getMessage());
    }

    private static List<String> names(List<String> supplied, int count, String prefix) {
        if (supplied == null) {
            List<String> result = new ArrayList<>(count);
            for (int index = 0; index < count; index++) result.add(prefix + (index + 1));
            return List.copyOf(result);
        }
        if (supplied.size() != count
                || supplied.stream().anyMatch(value -> value == null || value.isBlank()))
            throw new IllegalArgumentException(
                "one nonblank name is required per scan item");
        return List.copyOf(supplied);
    }

    private static double[] nan(int size) {
        double[] result = new double[size];
        Arrays.fill(result, Double.NaN);
        return result;
    }

    @FunctionalInterface private interface Chunk { void run(int chunk); }
    private record WeightedBase(double[] weightedDesign, double[] informationInverse,
            double[] responseResidual, double[] squareRootWeights,
            double residualSumSquares, int rank) { }
    private record PreparedPredictors(double[] values, boolean[] valid) { }
}
