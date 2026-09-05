/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.function.Supplier;
import org.jlinalg.loess.Loess;
import org.jlinalg.loess.LoessFamily;
import org.jlinalg.loess.LoessOptions;
import org.jlinalg.loess.LoessResult;

/** Deterministic benchmark shared with stats::loess direct surface. */
public final class LoessBenchmark {
    private static volatile double checksum;
    private LoessBenchmark() { }

    public static void main(String[] arguments) {
        int observations = property("jlinalg.benchmark.rows", 5_000);
        int warmups = property("jlinalg.benchmark.warmups", 2);
        int measurements = property("jlinalg.benchmark.measurements", 5);
        double span = doubleProperty("jlinalg.benchmark.span", 0.2);
        Data data = data(observations);
        LoessOptions options = new LoessOptions(
            span, 2, LoessFamily.GAUSSIAN, 4);
        Loess.Prepared prepared = Loess.prepare(data.x(), options);
        System.out.println("runtime,rows,span,degree,median_seconds,rows_per_second,checksum");
        measure("JLinAlg_end_to_end_direct", observations, span,
            warmups, measurements,
            () -> Loess.fit(data.x(), data.y(), null, options));
        measure("JLinAlg_prepared_direct", observations, span,
            warmups, measurements, () -> prepared.fit(data.y(), null));
    }

    private static void measure(String runtime, int observations, double span,
            int warmups, int measurements, Supplier<LoessResult> operation) {
        for (int iteration = 0; iteration < warmups; iteration++)
            consume(operation.get());
        double[] seconds = new double[measurements];
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime();
            consume(operation.get());
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
        }
        Arrays.sort(seconds);
        double median = seconds[seconds.length / 2];
        System.out.printf(java.util.Locale.ROOT,
            "%s,%d,%.6f,2,%.6f,%.2f,%.12f%n", runtime,
            observations, span, median, observations / median, checksum);
    }

    private static void consume(LoessResult result) {
        double[] fitted = result.fittedValues();
        checksum = fitted[0] + fitted[fitted.length / 2]
            + fitted[fitted.length - 1];
    }

    private static Data data(int observations) {
        double[] x = new double[observations];
        double[] y = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = -10.0 + 20.0 * row / (observations - 1.0);
            y[row] = Math.sin(0.7 * x[row]) + 0.05 * x[row]
                + 0.1 * Math.cos(0.37 * (row + 1.0));
        }
        return new Data(x, y);
    }

    private static int property(String name, int fallback) {
        return Integer.parseInt(System.getProperty(name,
            Integer.toString(fallback)));
    }

    private static double doubleProperty(String name, double fallback) {
        return Double.parseDouble(System.getProperty(name,
            Double.toString(fallback)));
    }

    private record Data(double[] x, double[] y) { }
}
