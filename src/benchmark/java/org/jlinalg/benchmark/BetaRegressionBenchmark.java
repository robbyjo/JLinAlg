/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.benchmark;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.BetaRegression;
import org.jlinalg.distributional.BetaRegressionResult;
import org.jlinalg.distributional.DistributionalFamilies;
import org.jlinalg.distributional.DistributionalModel;
import org.jlinalg.distributional.DistributionalOptions;
import org.jlinalg.distributional.DistributionalResult;
import org.jlinalg.gam.PenalizedPredictor;

/** Reproducible specialized-versus-generic beta-regression benchmark. */
public final class BetaRegressionBenchmark {
    private BetaRegressionBenchmark() { }

    public static void main(String[] arguments) {
        int rows = Integer.getInteger(
            "jlinalg.benchmark.beta.rows", 100_000);
        int repetitions = Integer.getInteger(
            "jlinalg.benchmark.beta.repetitions", 5);
        Data data = data(rows);
        PenalizedPredictor mean = PenalizedPredictor.linear(data.meanDesign());
        PenalizedPredictor precision = PenalizedPredictor.linear(
            data.precisionDesign());

        fitFast(data);
        fitGeneric(data, mean, precision);
        double[] fastTimes = new double[repetitions];
        double[] genericTimes = new double[repetitions];
        BetaRegressionResult fast = null;
        DistributionalResult generic = null;
        for (int repetition = 0; repetition < repetitions; repetition++) {
            long start = System.nanoTime();
            fast = fitFast(data);
            fastTimes[repetition] = elapsedMillis(start);
            start = System.nanoTime();
            generic = fitGeneric(data, mean, precision);
            genericTimes[repetition] = elapsedMillis(start);
        }
        if (fast == null || generic == null
                || !fast.converged() || !generic.converged()) {
            throw new IllegalStateException("benchmark fits must converge");
        }
        double fastMedian = median(fastTimes);
        double genericMedian = median(genericTimes);
        double coefficientDifference = maximumDifference(
            fast.meanCoefficients(),
            generic.parameter("mu").coefficients());
        double precisionDifference = maximumDifference(
            fast.precisionCoefficients(),
            generic.parameter("precision").coefficients());
        System.out.printf("rows=%d repetitions=%d%n", rows, repetitions);
        System.out.printf("specialized median_ms=%.3f iterations=%d%n",
            fastMedian, fast.iterations());
        System.out.printf("generic median_ms=%.3f iterations=%d%n",
            genericMedian, generic.iterations());
        System.out.printf("speedup=%.2fx%n", genericMedian / fastMedian);
        System.out.printf(
            "max_abs_mean_difference=%.3g max_abs_precision_difference=%.3g%n",
            coefficientDifference, precisionDifference);
        System.out.printf("coefficients=%s,%s%n",
            joined(fast.meanCoefficients()),
            joined(fast.precisionCoefficients()));
    }

    private static BetaRegressionResult fitFast(Data data) {
        return BetaRegression.fit(data.response(), data.meanDesign(),
            data.precisionDesign(),
            org.jlinalg.distributional.BetaRegressionOptions
                .variablePrecisionDefaults(),
            BackendPolicy.CPU);
    }

    private static DistributionalResult fitGeneric(
            Data data,
            PenalizedPredictor mean,
            PenalizedPredictor precision) {
        return DistributionalModel.fit(data.response(),
            List.of(mean, precision),
            DistributionalFamilies.betaMeanPrecision(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
    }

    private static Data data(int rows) {
        double[] response = new double[rows];
        double[][] mean = new double[rows][4];
        double[][] precision = new double[rows][2];
        for (int row = 0; row < rows; row++) {
            double x1 = -1.0 + 2.0 * row / Math.max(1.0, rows - 1.0);
            double x2 = Math.sin(0.017 * row);
            double x3 = Math.cos(0.031 * row);
            mean[row][0] = 1.0;
            mean[row][1] = x1;
            mean[row][2] = x2;
            mean[row][3] = x3;
            precision[row][0] = 1.0;
            precision[row][1] = x1;
            double mu = logistic(-0.4 + 0.8 * x1 - 0.35 * x2 + 0.2 * x3);
            double phi = Math.exp(3.8 + 0.25 * x1);
            double standardized = (Math.sin(1.73 * row + 0.2)
                + 0.55 * Math.cos(0.47 * row)) / 1.14;
            double standardDeviation = Math.sqrt(mu * (1.0 - mu)
                / (phi + 1.0));
            response[row] = Math.max(1e-5,
                Math.min(1.0 - 1e-5,
                    mu + standardDeviation * standardized));
        }
        return new Data(response, mean, precision);
    }

    private static double logistic(double value) {
        return value >= 0.0 ? 1.0 / (1.0 + Math.exp(-value))
            : Math.exp(value) / (1.0 + Math.exp(value));
    }

    private static double elapsedMillis(long start) {
        return (System.nanoTime() - start) / 1_000_000.0;
    }

    private static double median(double[] values) {
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        return sorted[sorted.length / 2];
    }

    private static double maximumDifference(double[] first, double[] second) {
        double maximum = 0.0;
        for (int index = 0; index < first.length; index++) {
            maximum = Math.max(maximum,
                Math.abs(first[index] - second[index]));
        }
        return maximum;
    }

    private static String joined(double[] values) {
        StringBuilder result = new StringBuilder();
        for (int index = 0; index < values.length; index++) {
            if (index > 0) result.append(',');
            result.append(String.format("%.16g", values[index]));
        }
        return result.toString();
    }

    private record Data(
            double[] response,
            double[][] meanDesign,
            double[][] precisionDesign) { }
}
