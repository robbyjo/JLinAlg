/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.Random;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.survival.CoxOptions;
import org.jlinalg.survival.CoxRegression;
import org.jlinalg.survival.CoxResult;
import org.jlinalg.survival.CoxScoreVariance;
import org.jlinalg.survival.CoxSurvivalData;
import org.jlinalg.survival.FastCoxAssociation;

/** Reproducible throughput exercise for the prepared survival paths. */
public final class CoxPipelineBenchmark {
    private static volatile double checksum;
    private CoxPipelineBenchmark() { }

    public static void main(String[] arguments) {
        int rows = property("jlinalg.benchmark.rows", 2_000);
        int variables = property("jlinalg.benchmark.variables", 512);
        int warmups = property("jlinalg.benchmark.warmups", 1);
        int measurements = property("jlinalg.benchmark.measurements", 3);
        int parallelism = property("jlinalg.benchmark.parallelism",
            Math.max(1, Runtime.getRuntime().availableProcessors()));
        Data data = data(rows, variables);
        AssociationEngineOptions execution = AssociationEngineOptions
            .cpuParallel().withParallelism(parallelism).withChunkSize(32);
        FastCoxAssociation prepared = FastCoxAssociation.prepare(
            data.right(), data.covariates(), null, CoxOptions.defaults(),
            execution);

        System.out.println("benchmark,rows,variables,parallelism,median_seconds,"
            + "variables_per_second");
        measure("prepared_cox_score", rows, variables, parallelism,
            warmups, measurements, () -> {
                AssociationBatchResult result = prepared.scan(data.predictors(),
                    null, execution, CoxScoreVariance.MODEL_BASED);
                checksum += result.beta()[variables - 1];
            }, variables);

        int exact = Math.min(32, variables);
        try (CoxRegression.Prepared fitter = CoxRegression.prepare(data.right(),
                CoxOptions.defaults(), BackendPolicy.CPU)) {
            measure("full_cox_refit", rows, exact, 1, warmups, measurements,
                () -> {
                    for (int variable = 0; variable < exact; variable++) {
                        double[][] design = append(data.covariates(),
                            data.predictors(), variable);
                        CoxResult result = fitter.fit(design, null);
                        checksum += result.beta()[design[0].length - 1];
                    }
                }, exact);
        }

        measure("start_stop_sweep_fit", rows, 1, 1, warmups, measurements,
            () -> {
                CoxResult result = CoxRegression.fit(data.counting(),
                    data.covariates(), null, CoxOptions.defaults(),
                    BackendPolicy.CPU);
                checksum += result.beta()[0];
            }, 1);
    }

    private static void measure(String name, int rows, int variables,
            int parallelism, int warmups, int measurements, Runnable operation,
            int completed) {
        for (int index = 0; index < warmups; index++) operation.run();
        double[] seconds = new double[measurements];
        for (int index = 0; index < measurements; index++) {
            long started = System.nanoTime();
            operation.run();
            seconds[index] = (System.nanoTime() - started) / 1e9;
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        System.out.printf(java.util.Locale.ROOT, "%s,%d,%d,%d,%.6f,%.2f%n",
            name, rows, variables, parallelism, median, completed / median);
    }

    private static Data data(int rows, int variables) {
        Random random = new Random(20260905L);
        double[][] covariates = new double[rows][2];
        double[][] predictors = new double[rows][variables];
        double[] stop = new double[rows];
        double[] start = new double[rows];
        boolean[] event = new boolean[rows];
        for (int row = 0; row < rows; row++) {
            covariates[row][0] = random.nextGaussian();
            covariates[row][1] = random.nextGaussian();
            stop[row] = 0.25 + 10.0 * random.nextDouble();
            start[row] = row % 3 == 0 ? stop[row] * random.nextDouble() * 0.7 : 0.0;
            event[row] = random.nextDouble() < 0.55;
            for (int variable = 0; variable < variables; variable++)
                predictors[row][variable] = random.nextGaussian();
        }
        event[0] = true;
        return new Data(CoxSurvivalData.rightCensored(stop, event),
            new CoxSurvivalData(start, stop, event, null), covariates,
            predictors);
    }

    private static double[][] append(double[][] covariates,
            double[][] predictors, int variable) {
        double[][] result = new double[covariates.length][3];
        for (int row = 0; row < result.length; row++) {
            result[row][0] = covariates[row][0];
            result[row][1] = covariates[row][1];
            result[row][2] = predictors[row][variable];
        }
        return result;
    }

    private static int property(String name, int fallback) {
        return Integer.parseInt(System.getProperty(name,
            Integer.toString(fallback)));
    }

    private record Data(CoxSurvivalData right, CoxSurvivalData counting,
        double[][] covariates, double[][] predictors) { }
}
