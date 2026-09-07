/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.Arrays;
import org.jlinalg.regression.KernelRegression;
import org.jlinalg.regression.MultinomialRegression;
import org.jlinalg.regression.MultivariateRegression;
import org.jlinalg.regression.PartiallyLinearRegression;
import org.jlinalg.regression.QuantileRegression;
import org.jlinalg.regression.SuperSmoother;

/** Reproducible JVM throughput profile for the remaining regression families. */
public final class RegressionFamiliesBenchmark {
    private static volatile double checksum;
    private RegressionFamiliesBenchmark() { }

    public static void main(String[] arguments) {
        int rows = property("jlinalg.benchmark.regression.rows", 512);
        int warmups = property("jlinalg.benchmark.regression.warmups", 1);
        int measurements = property("jlinalg.benchmark.regression.measurements", 3);
        Data data = data(rows);
        System.out.println("operation,rows,median_seconds,rows_per_second,checksum");
        measure("multivariate_ols", rows, warmups, measurements,
            () -> MultivariateRegression.fit(data.multivariateResponse(), data.design()));
        measure("multinomial_logit", rows, warmups, measurements,
            () -> MultinomialRegression.fit(data.classes(), data.design(), 3));
        measure("quantile", rows, warmups, measurements,
            () -> QuantileRegression.fit(data.response(), data.design(), 0.5));
        measure("kernel", rows, warmups, measurements,
            () -> KernelRegression.fit(data.predictor(), data.response()));
        measure("supersmoother", rows, warmups, measurements,
            () -> SuperSmoother.fit(data.predictor(), data.response()));
        measure("partially_linear", rows, warmups, measurements,
            () -> PartiallyLinearRegression.fit(data.response(), data.design(), data.predictor(), 0.4));
    }

    private static void measure(String name, int rows, int warmups, int measurements, Runnable operation) {
        for (int iteration = 0; iteration < warmups; iteration++) operation.run();
        double[] seconds = new double[measurements];
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime(); operation.run();
            seconds[iteration] = (System.nanoTime() - started) / 1e9;
        }
        Arrays.sort(seconds); double median = seconds[seconds.length / 2];
        System.out.printf(java.util.Locale.ROOT, "%s,%d,%.6f,%.2f,%.12f%n",
            name, rows, median, rows / median, checksum);
    }

    private static Data data(int rows) {
        double[] predictor = new double[rows], response = new double[rows];
        double[][] design = new double[rows][2], multivariateResponse = new double[rows][2];
        int[] classes = new int[rows];
        for (int row = 0; row < rows; row++) {
            predictor[row] = -3.0 + 6.0 * row / (rows - 1.0);
            design[row][0] = 1.0; design[row][1] = predictor[row];
            response[row] = 1.0 + 0.4 * predictor[row] + Math.sin(predictor[row]);
            multivariateResponse[row][0] = response[row];
            multivariateResponse[row][1] = -0.5 + 0.8 * predictor[row];
            classes[row] = predictor[row] < -1 ? 0 : predictor[row] > 1 ? 2 : 1;
        }
        return new Data(predictor, response, design, multivariateResponse, classes);
    }

    private static int property(String name, int fallback) {
        return Integer.parseInt(System.getProperty(name, Integer.toString(fallback)));
    }

    private record Data(double[] predictor, double[] response, double[][] design,
                        double[][] multivariateResponse, int[] classes) { }
}
