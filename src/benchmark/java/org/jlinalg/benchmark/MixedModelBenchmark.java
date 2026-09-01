/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.LinearMixedModel;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.reml.RemlOptions;

/** Dense-versus-sparse random-intercept end-to-end macrobenchmark. */
public final class MixedModelBenchmark {
    private static volatile double checksum;

    private MixedModelBenchmark() { }

    public static void main(String[] arguments) {
        Logger.getLogger("jdistlib.accelerator.ComputeBackends")
            .setLevel(Level.WARNING);
        int rows = integerProperty("jlinalg.benchmark.rows", 10_000);
        int groups = integerProperty("jlinalg.benchmark.groups",
            Math.max(2, rows / 5));
        int repeats = integerProperty("jlinalg.benchmark.measurements", 3);
        Data data = data(rows, groups);
        RemlOptions options = RemlOptions.builder()
            .initialVariances(1.0, 1.0).maximumIterations(30).build();
        System.out.println("benchmark,rows,groups,seconds,rows_per_second,"
            + "equation_nnz,factor_nnz");
        double[] times = new double[repeats];
        int equationNonzeros = 0;
        int factorNonzeros = 0;
        for (int iteration = 0; iteration < repeats; iteration++) {
            long started = System.nanoTime();
            var fit = SparseLinearMixedModel.fit(data.response(), data.fixed(),
                List.of(data.groups()), options, BackendPolicy.CPU);
            times[iteration] = (System.nanoTime() - started) / 1e9;
            checksum += fit.beta()[0];
            equationNonzeros = fit.equationNonzeroCount();
            factorNonzeros = fit.factorNonzeroCount();
        }
        java.util.Arrays.sort(times);
        double median = times[times.length / 2];
        System.out.printf(java.util.Locale.ROOT,
            "sparse_reml,%d,%d,%.6f,%.2f,%d,%d%n", rows, groups,
            median, rows / median, equationNonzeros, factorNonzeros);

        if (rows <= 2_000) {
            long started = System.nanoTime();
            var dense = LinearMixedModel.fit(data.response(), data.fixed(),
                List.of(data.groups()), options, BackendPolicy.CPU);
            double seconds = (System.nanoTime() - started) / 1e9;
            checksum += dense.beta()[0];
            System.out.printf(java.util.Locale.ROOT,
                "dense_reml,%d,%d,%.6f,%.2f,NA,NA%n", rows, groups,
                seconds, rows / seconds);
        }
    }

    private static Data data(int rows, int groupCount) {
        Random random = new Random(20260901L);
        double[] groupEffects = new double[groupCount];
        for (int group = 0; group < groupCount; group++)
            groupEffects[group] = random.nextGaussian();
        List<String> groups = new ArrayList<>(rows);
        double[] response = new double[rows];
        double[][] fixed = new double[rows][2];
        for (int row = 0; row < rows; row++) {
            int group = row % groupCount;
            double x = random.nextGaussian();
            groups.add("g" + group);
            fixed[row][0] = 1.0;
            fixed[row][1] = x;
            response[row] = 2.0 + 0.4 * x + groupEffects[group]
                + 0.7 * random.nextGaussian();
        }
        return new Data(response, fixed,
            RandomEffectTerm.randomIntercept("group", groups));
    }

    private static int integerProperty(String name, int defaultValue) {
        int value = Integer.getInteger(name, defaultValue);
        if (value < 1) throw new IllegalArgumentException(name + " must be positive");
        return value;
    }

    private record Data(double[] response, double[][] fixed,
                        RandomEffectTerm groups) { }
}
