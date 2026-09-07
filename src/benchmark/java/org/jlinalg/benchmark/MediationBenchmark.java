/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mediation.MediationMixedResult;
import org.jlinalg.mediation.MediationResult;
import org.jlinalg.mediation.MediationAnalysis;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.reml.RemlOptions;

/** Deterministic Gaussian mediation benchmark shared with the R comparator. */
public final class MediationBenchmark {
    private static volatile double checksum;

    private MediationBenchmark() { }

    public static void main(String[] arguments) {
        int rows = integerProperty("jlinalg.benchmark.mediation.rows", 2_000);
        int groups = integerProperty("jlinalg.benchmark.mediation.groups", 100);
        int warmups = integerProperty("jlinalg.benchmark.mediation.warmups", 1);
        int measurements = integerProperty(
            "jlinalg.benchmark.mediation.measurements", 3);
        Data data = data(rows, groups);
        System.out.println("runtime,model,rows,groups,median_seconds,"
            + "rows_per_second,a,b,indirect,direct,total,checksum");

        measure("JLinAlg", "ols", rows, groups, warmups, measurements,
            () -> MediationAnalysis.fit(
                data.outcome(), data.treatment(), data.mediator(),
                data.covariates(),
                org.jlinalg.ols.OlsOptions.defaults(), BackendPolicy.CPU));

        RemlOptions options = RemlOptions.builder()
            .initialVariances(0.5, 0.5)
            .maximumIterations(40)
            .build();
        measure("JLinAlg", "sparse_mixed", rows, groups, warmups,
            measurements, () -> MediationAnalysis.fitMixed(
                data.outcome(), data.treatment(), data.mediator(),
                data.covariates(), List.of(data.randomIntercept()), options,
                BackendPolicy.CPU));
    }

    private static void measure(
            String runtime,
            String model,
            int rows,
            int groups,
            int warmups,
            int measurements,
            Supplier<Object> operation) {
        for (int iteration = 0; iteration < warmups; iteration++)
            consume(operation.get());
        double[] seconds = new double[measurements];
        Object result = null;
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime();
            result = operation.get();
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
            consume(result);
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        MediationEffectValues effects = effects(result);
        System.out.printf(java.util.Locale.ROOT,
            "%s,%s,%d,%d,%.6f,%.2f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f%n",
            runtime, model, rows, groups, median, rows / median,
            effects.a(), effects.b(), effects.indirect(), effects.direct(),
            effects.total(), checksum);
    }

    private static void consume(Object result) {
        MediationEffectValues effects = effects(result);
        checksum = effects.a() + effects.b() + effects.indirect()
            + effects.direct() + effects.total();
    }

    private static MediationEffectValues effects(Object result) {
        if (result instanceof MediationResult mediation) {
            return new MediationEffectValues(
                mediation.aPath().estimate(), mediation.bPath().estimate(),
                mediation.indirectEffect().estimate(),
                mediation.directEffect().estimate(),
                mediation.totalEffect().estimate());
        }
        MediationMixedResult mediation = (MediationMixedResult) result;
        return new MediationEffectValues(
            mediation.aPath().estimate(), mediation.bPath().estimate(),
            mediation.indirectEffect().estimate(),
            mediation.directEffect().estimate(),
            mediation.totalEffect().estimate());
    }

    private static Data data(int rows, int groupCount) {
        double[] treatment = new double[rows];
        double[] mediator = new double[rows];
        double[] outcome = new double[rows];
        double[][] covariates = new double[rows][1];
        List<String> groups = new ArrayList<>(rows);
        for (int row = 0; row < rows; row++) {
            double index = row + 1.0;
            double x = -1.0 + 2.0 * row / (rows - 1.0);
            double covariate = Math.sin(0.17 * index);
            double groupEffect = 0.4 * Math.sin(0.7 * (row % groupCount));
            treatment[row] = x;
            covariates[row][0] = covariate;
            mediator[row] = 0.8 + 0.6 * x + 0.4 * covariate
                + groupEffect + 0.1 * Math.cos(0.31 * index);
            outcome[row] = 0.3 + 0.25 * x + 1.1 * mediator[row]
                + 0.2 * covariate + groupEffect
                + 0.1 * Math.sin(0.23 * index);
            groups.add("g" + (row % groupCount));
        }
        return new Data(outcome, treatment, mediator, covariates,
            RandomEffectTerm.randomIntercept("group", groups));
    }

    private static int integerProperty(String name, int fallback) {
        int value = Integer.getInteger(name, fallback);
        if (value < 1) throw new IllegalArgumentException(
            name + " must be positive");
        return value;
    }

    private record Data(
            double[] outcome,
            double[] treatment,
            double[] mediator,
            double[][] covariates,
            RandomEffectTerm randomIntercept) { }

    private record MediationEffectValues(
            double a, double b, double indirect, double direct, double total) { }
}
