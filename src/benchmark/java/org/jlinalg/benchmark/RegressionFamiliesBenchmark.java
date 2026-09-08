/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.function.DoubleSupplier;
import org.jlinalg.compute.BackendPolicy;
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
        if(rows<16||warmups<0||measurements<1)throw new IllegalArgumentException("invalid benchmark controls");
        Data data = data(rows);
        System.out.println("operation,rows,median_seconds,rows_per_second,checksum");
        measure("multivariate_ols", rows, warmups, measurements,
            () -> Arrays.stream(MultivariateRegression.fit(data.multivariateResponse(), data.design(),BackendPolicy.CPU).coefficients()).sum());
        measure("multinomial_logit", rows, warmups, measurements,
            () -> {var fit=MultinomialRegression.fit(data.classes(),data.design(),3);if(!fit.converged())throw new IllegalStateException("multinomial failed accuracy gate");return fit.logLikelihood();});
        measure("quantile", rows, warmups, measurements,
            () -> {var fit=QuantileRegression.fit(data.response(),data.design(),.5);if(!fit.converged())throw new IllegalStateException("quantile failed accuracy gate");return Arrays.stream(fit.coefficients()).sum();});
        measure("kernel", rows, warmups, measurements,
            () -> Arrays.stream(KernelRegression.fit(data.predictor(),data.response(),.4).fittedValues()).sum());
        measure("supersmoother", rows, warmups, measurements,
            () -> Arrays.stream(SuperSmoother.fit(data.predictor(),data.response()).fittedValues()).sum());
        measure("partially_linear", rows, warmups, measurements,
            () -> Arrays.stream(PartiallyLinearRegression.fitWithInference(data.response(),data.design(),data.predictor(),.4,BackendPolicy.CPU).fit().coefficients()).sum());
    }

    private static void measure(String name, int rows, int warmups, int measurements, DoubleSupplier operation) {
        for (int iteration = 0; iteration < warmups; iteration++) checksum=operation.getAsDouble();
        double[] seconds = new double[measurements];
        for (int iteration = 0; iteration < measurements; iteration++) {
            long started = System.nanoTime(); checksum=operation.getAsDouble();
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
            double linear=Math.cos(row*1.7)+.2*predictor[row];
            design[row][0] = 1.0; design[row][1] = linear;
            response[row] = 1.0 + 0.4 * linear + Math.sin(predictor[row])+.15*Math.cos(row*.37);
            multivariateResponse[row][0] = response[row];
            multivariateResponse[row][1] = -0.5 + 0.8 * linear+.1*Math.sin(row*.8);
            double l1=Math.exp(.3+.5*linear),l2=Math.exp(-.2-.3*linear),u=((row*37+11)%101+.5)/101;
            classes[row]=u<1/(1+l1+l2)?0:u<(1+l1)/(1+l1+l2)?1:2;
        }
        return new Data(predictor, response, design, multivariateResponse, classes);
    }

    private static int property(String name, int fallback) {
        return Integer.parseInt(System.getProperty(name, Integer.toString(fallback)));
    }

    private record Data(double[] predictor, double[] response, double[][] design,
                        double[][] multivariateResponse, int[] classes) { }
}
