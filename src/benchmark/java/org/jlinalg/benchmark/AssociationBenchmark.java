/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.association.AssociationModels;
import org.jlinalg.association.FastGlmAssociation;
import org.jlinalg.association.FastOlsAssociation;
import org.jlinalg.association.ParallelAssociationEngine;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.ols.OlsOptions;

/** Deterministic end-to-end association throughput benchmark. */
public final class AssociationBenchmark {
    private static volatile double checksum;

    private AssociationBenchmark() { }

    public static void main(String[] arguments) {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        int observations = integerProperty("jlinalg.benchmark.rows", 10_000);
        int variables = integerProperty("jlinalg.benchmark.variables", 2_000);
        int covariates = integerProperty("jlinalg.benchmark.covariates", 6);
        int warmups = integerProperty("jlinalg.benchmark.warmups", 2);
        int measurements = integerProperty("jlinalg.benchmark.measurements", 5);
        int parallelism = integerProperty("jlinalg.benchmark.parallelism",
            Math.max(1, Runtime.getRuntime().availableProcessors()));
        Data data = data(observations, variables, covariates);
        AssociationEngineOptions execution =
            AssociationEngineOptions.cpuParallel()
                .withParallelism(parallelism).withChunkSize(128);
        FastGlmAssociation logistic = FastGlmAssociation.prepare(
            data.binary(), data.covariates(), GlmFamilies.binomial(),
            null, null, GlmOptions.defaults(), execution);

        System.out.println("benchmark,rows,variables,covariates,parallelism,"
            + "median_seconds,variables_per_second");
        measure("fast_ols", observations, variables, covariates,
            parallelism, warmups, measurements,
            () -> FastOlsAssociation.scanPredictors(data.gaussian(),
                data.covariates(), data.predictors(), null, null, null,
                OlsOptions.defaults(), execution), variables);
        measure("prepared_glm_score", observations, variables, covariates,
            parallelism, warmups, measurements,
            () -> logistic.scan(data.predictors(), null, execution), variables);

        int exactVariables = Math.min(variables, 100);
        double[][] exactPredictors = firstColumns(
            data.predictors(), exactVariables);
        measure("exact_parallel_ols", observations, exactVariables, covariates,
            parallelism, warmups, measurements,
            () -> ParallelAssociationEngine.scanPredictors(data.gaussian(),
                data.covariates(), exactPredictors, null,
                AssociationModels.ols(OlsOptions.defaults()), execution),
            exactVariables);
    }

    private static void measure(
            String name, int rows, int variables, int covariates,
            int parallelism, int warmups, int measurements,
            Operation operation, int completedVariables) {
        for (int iteration = 0; iteration < warmups; iteration++)
            consume(operation.run());
        double[] seconds = new double[measurements];
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime();
            consume(operation.run());
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        System.out.printf(java.util.Locale.ROOT,
            "%s,%d,%d,%d,%d,%.6f,%.2f%n", name, rows, variables,
            covariates, parallelism, median, completedVariables / median);
    }

    private static void consume(AssociationBatchResult result) {
        double[] beta = result.beta();
        checksum += beta.length == 0 ? 0.0 : beta[beta.length - 1];
    }

    private static Data data(int rows, int variables, int columns) {
        if (rows <= columns + 2 || variables < 1 || columns < 1)
            throw new IllegalArgumentException("benchmark dimensions are invalid");
        Random random = new Random(20260901L);
        double[][] covariates = new double[rows][columns];
        double[][] predictors = new double[rows][variables];
        double[] gaussian = new double[rows];
        double[] binary = new double[rows];
        for (int row = 0; row < rows; row++) {
            covariates[row][0] = 1.0;
            for (int column = 1; column < columns; column++)
                covariates[row][column] = random.nextGaussian();
            double eta = -0.3 + 0.15 * covariates[row][Math.min(1, columns - 1)];
            gaussian[row] = 1.5 + 0.4 * eta + random.nextGaussian();
            double probability = 1.0 / (1.0 + Math.exp(-eta));
            binary[row] = random.nextDouble() < probability ? 1.0 : 0.0;
            for (int variable = 0; variable < variables; variable++)
                predictors[row][variable] = dosage(random);
        }
        return new Data(gaussian, binary, covariates, predictors);
    }

    private static double dosage(Random random) {
        return (random.nextDouble() < 0.28 ? 1.0 : 0.0)
            + (random.nextDouble() < 0.28 ? 1.0 : 0.0);
    }

    private static double[][] firstColumns(double[][] matrix, int count) {
        double[][] result = new double[matrix.length][count];
        for (int row = 0; row < matrix.length; row++)
            System.arraycopy(matrix[row], 0, result[row], 0, count);
        return result;
    }

    private static int integerProperty(String name, int defaultValue) {
        int value = Integer.getInteger(name, defaultValue);
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }

    @FunctionalInterface
    private interface Operation { AssociationBatchResult run(); }
    private record Data(double[] gaussian, double[] binary,
                        double[][] covariates, double[][] predictors) { }
}
