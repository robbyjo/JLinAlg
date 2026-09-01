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
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;

/** Ordered parallel execution for changing predictors or responses. */
public final class ParallelAssociationEngine {
    private ParallelAssociationEngine() { }

    /**
     * Appends each observation-by-variable predictor column to a fixed design
     * and extracts association statistics for that final column.
     */
    public static AssociationBatchResult scanPredictors(
            double[] response,
            double[][] baseDesign,
            double[][] predictors,
            List<String> predictorNames,
            AssociationFitter fitter,
            AssociationEngineOptions options) {
        if (response == null) throw new IllegalArgumentException("response is required");
        double[] base = MatrixOps.rowMajor(baseDesign, response.length);
        double[] variables = MatrixOps.rowMajorUnchecked(predictors, response.length);
        return scanPredictors(response, base, response.length, baseDesign[0].length,
            variables, predictors[0].length, predictorNames, fitter, options);
    }

    /** Contiguous row-major predictor scan. */
    public static AssociationBatchResult scanPredictors(
            double[] response,
            double[] baseDesign,
            int rows,
            int baseColumns,
            double[] predictors,
            int predictorCount,
            List<String> predictorNames,
            AssociationFitter fitter,
            AssociationEngineOptions options) {
        MatrixOps.validateModelData(response, baseDesign, rows, baseColumns);
        if (predictors == null || predictorCount < 1
                || predictors.length != rows * predictorCount
                || fitter == null || options == null) {
            throw new IllegalArgumentException("predictor data, fitter, or options are invalid");
        }
        List<String> names = names(predictorNames, predictorCount, "variable");
        Accumulator result = new Accumulator(predictorCount);
        long started = System.nanoTime();
        int chunks = (predictorCount + options.chunkSize() - 1) / options.chunkSize();
        execute(chunks, options.parallelism(), chunk -> {
            double[] augmented = new double[rows * (baseColumns + 1)];
            for (int variable = chunk * options.chunkSize();
                    variable < Math.min(predictorCount,
                        (chunk + 1) * options.chunkSize()); variable++) {
                try {
                    appendPredictor(baseDesign, baseColumns, predictors,
                        predictorCount, variable, rows, augmented,
                        options.predictorMissingPolicy());
                    AssociationStatistics statistics = fitter.fit(
                        response, augmented, rows, baseColumns + 1,
                        options.backendPolicy());
                    result.set(variable, statistics, baseColumns);
                } catch (RuntimeException exception) {
                    handleFailure(result, variable, names.get(variable),
                        exception, options.failurePolicy());
                }
            }
        });
        return result.toResult(names, baseColumns,
            Math.min(options.parallelism(), chunks),
            System.nanoTime() - started);
    }

    /**
     * Fits each observation-by-trait response column against one fixed design
     * and extracts the caller-selected fixed-effect coefficient.
     */
    public static AssociationBatchResult scanResponses(
            double[][] responses,
            double[][] design,
            int coefficientIndex,
            List<String> responseNames,
            AssociationFitter fitter,
            AssociationEngineOptions options) {
        if (responses == null || responses.length == 0)
            throw new IllegalArgumentException("responses are required");
        double[] responseMatrix = MatrixOps.rowMajorUnchecked(responses, responses.length);
        double[] fixed = MatrixOps.rowMajor(design, responses.length);
        return scanResponses(responseMatrix, responses.length, responses[0].length,
            fixed, design[0].length, coefficientIndex,
            responseNames, fitter, options);
    }

    /** Contiguous row-major response scan. */
    public static AssociationBatchResult scanResponses(
            double[] responses,
            int rows,
            int responseCount,
            double[] design,
            int columns,
            int coefficientIndex,
            List<String> responseNames,
            AssociationFitter fitter,
            AssociationEngineOptions options) {
        if (responses == null || responseCount < 1
                || responses.length != rows * responseCount
                || coefficientIndex < 0 || coefficientIndex >= columns
                || fitter == null || options == null) {
            throw new IllegalArgumentException("response scan dimensions or controls are invalid");
        }
        MatrixOps.validateModelData(new double[rows], design, rows, columns);
        List<String> names = names(responseNames, responseCount, "response");
        Accumulator result = new Accumulator(responseCount);
        long started = System.nanoTime();
        int chunks = (responseCount + options.chunkSize() - 1) / options.chunkSize();
        execute(chunks, options.parallelism(), chunk -> {
            double[] response = new double[rows];
            for (int trait = chunk * options.chunkSize();
                    trait < Math.min(responseCount,
                        (chunk + 1) * options.chunkSize()); trait++) {
                for (int row = 0; row < rows; row++)
                    response[row] = responses[row * responseCount + trait];
                try {
                    AssociationStatistics statistics = fitter.fit(
                        response, design, rows, columns, options.backendPolicy());
                    result.set(trait, statistics, coefficientIndex);
                } catch (RuntimeException exception) {
                    handleFailure(result, trait, names.get(trait),
                        exception, options.failurePolicy());
                }
            }
        });
        return result.toResult(names, coefficientIndex,
            Math.min(options.parallelism(), chunks),
            System.nanoTime() - started);
    }

    private static void appendPredictor(
            double[] base, int baseColumns,
            double[] predictors, int predictorCount, int predictor,
            int rows, double[] augmented, VariableMissingPolicy missingPolicy) {
        double sum = 0.0;
        int finite = 0;
        for (int row = 0; row < rows; row++) {
            double value = predictors[row * predictorCount + predictor];
            if (Double.isFinite(value)) { sum += value; finite++; }
            else if (missingPolicy == VariableMissingPolicy.ERROR)
                throw new IllegalArgumentException(
                    "non-finite predictor at row " + row);
        }
        if (finite == 0) throw new IllegalArgumentException("predictor has no finite values");
        double mean = sum / finite;
        int columns = baseColumns + 1;
        for (int row = 0; row < rows; row++) {
            System.arraycopy(base, row * baseColumns,
                augmented, row * columns, baseColumns);
            double value = predictors[row * predictorCount + predictor];
            augmented[row * columns + baseColumns] =
                Double.isFinite(value) ? value : mean;
        }
    }

    private static void execute(
            int chunks, int parallelism, ChunkOperation operation) {
        ForkJoinPool pool = new ForkJoinPool(Math.min(parallelism, chunks));
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
        } finally {
            pool.shutdownNow();
        }
    }

    private static void handleFailure(
            Accumulator result, int index, String name,
            RuntimeException exception, AssociationFailurePolicy policy) {
        if (policy == AssociationFailurePolicy.FAIL_FAST) throw exception;
        result.failures.add(new AssociationFailure(index, name,
            exception.getClass().getSimpleName(),
            exception.getMessage() == null ? "" : exception.getMessage()));
    }

    private static List<String> names(
            List<String> supplied, int count, String prefix) {
        if (supplied == null) {
            List<String> result = new ArrayList<>(count);
            for (int index = 0; index < count; index++) result.add(prefix + (index + 1));
            return List.copyOf(result);
        }
        if (supplied.size() != count
                || supplied.stream().anyMatch(value -> value == null || value.isBlank()))
            throw new IllegalArgumentException("one nonblank name is required per scan item");
        return List.copyOf(supplied);
    }

    @FunctionalInterface
    private interface ChunkOperation { void run(int chunk); }

    private static final class Accumulator {
        private final double[] beta;
        private final double[] standardErrors;
        private final double[] statistics;
        private final double[] pValues;
        private final double[] negativeLog10PValues;
        private final double[] degreesOfFreedom;
        private final ConcurrentLinkedQueue<AssociationFailure> failures =
            new ConcurrentLinkedQueue<>();

        Accumulator(int size) {
            beta = nan(size);
            standardErrors = nan(size);
            statistics = nan(size);
            pValues = nan(size);
            negativeLog10PValues = nan(size);
            degreesOfFreedom = nan(size);
        }

        void set(int destination, AssociationStatistics source, int coefficient) {
            double[] sourceBeta = source.beta();
            if (coefficient >= sourceBeta.length)
                throw new IllegalArgumentException(
                    "requested coefficient is absent from fitter result");
            beta[destination] = sourceBeta[coefficient];
            standardErrors[destination] = source.standardErrors()[coefficient];
            statistics[destination] = source.statistics()[coefficient];
            pValues[destination] = source.pValues()[coefficient];
            negativeLog10PValues[destination] =
                source.negativeLog10PValues()[coefficient];
            degreesOfFreedom[destination] = source.degreesOfFreedom()[coefficient];
        }

        AssociationBatchResult toResult(
                List<String> names, int coefficient,
                int parallelism, long elapsed) {
            List<AssociationFailure> ordered = failures.stream()
                .sorted(Comparator.comparingInt(AssociationFailure::index)).toList();
            return new AssociationBatchResult(names, beta, standardErrors,
                statistics, pValues, degreesOfFreedom, negativeLog10PValues,
                ordered,
                coefficient, parallelism, elapsed);
        }

        private static double[] nan(int size) {
            double[] result = new double[size];
            Arrays.fill(result, Double.NaN);
            return result;
        }
    }
}
