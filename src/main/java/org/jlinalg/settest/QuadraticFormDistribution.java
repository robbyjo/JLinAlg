/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.Arrays;
import jdistlib.ChiSquare;

/** Deterministic positive chi-square-mixture survival probabilities. */
final class QuadraticFormDistribution {
    private static final int INTERVALS_PER_BLOCK = 128;
    private static final int MAXIMUM_BLOCKS = 512;
    private static final int REQUIRED_QUIET_BLOCKS = 6;
    private static final double ABSOLUTE_TOLERANCE = 1e-9;

    private QuadraticFormDistribution() { }

    record Tail(double pValue, String method,
        double degreesOfFreedom, double scale) { }

    static Tail survival(double statistic, double[] suppliedEigenvalues) {
        if (!Double.isFinite(statistic) || statistic < 0
                || suppliedEigenvalues == null)
            throw new IllegalArgumentException(
                "quadratic-form statistic and eigenvalues are invalid");
        double[] eigenvalues = Arrays.stream(suppliedEigenvalues)
            .filter(value -> value > 1e-12 && Double.isFinite(value))
            .toArray();
        if (eigenvalues.length == 0)
            throw new IllegalArgumentException(
                "quadratic-form kernel has no positive eigenvalues");
        double sum = Arrays.stream(eigenvalues).sum();
        double sumSquares = Arrays.stream(eigenvalues)
            .map(value -> value * value).sum();
        double degrees = sum * sum / sumSquares;
        double scale = sumSquares / sum;
        if (eigenvalues.length == 1)
            return new Tail(probability(ChiSquare.cumulative(
                statistic / eigenvalues[0], 1, false, false)),
                "exact-scaled-chi-square", 1, eigenvalues[0]);
        boolean equal = true;
        for (int index = 1; index < eigenvalues.length; index++)
            if (Math.abs(eigenvalues[index] - eigenvalues[0])
                    > 1e-12 * Math.max(eigenvalues[index], eigenvalues[0])) {
                equal = false;
                break;
            }
        if (equal)
            return new Tail(probability(ChiSquare.cumulative(
                statistic / eigenvalues[0], eigenvalues.length,
                false, false)), "exact-equal-eigenvalue-chi-square",
                eigenvalues.length, eigenvalues[0]);
        Double imhof = imhof(statistic, eigenvalues);
        if (imhof != null)
            return new Tail(probability(imhof), "imhof", degrees, scale);
        return new Tail(probability(ChiSquare.cumulative(
            statistic / scale, degrees, false, false)),
            "satterthwaite-fallback", degrees, scale);
    }

    static double critical(double[] eigenvalues, double survivalProbability) {
        if (!(survivalProbability > 0 && survivalProbability <= 1)
                || !Double.isFinite(survivalProbability))
            throw new IllegalArgumentException(
                "survival probability must be in (0,1]");
        if (survivalProbability == 1) return 0;
        double high = Arrays.stream(eigenvalues).sum();
        for (int iteration = 0; iteration < 128
                && survival(high, eigenvalues).pValue() > survivalProbability;
                iteration++)
            high *= 2;
        double low = 0;
        for (int iteration = 0; iteration < 64; iteration++) {
            double middle = (low + high) * 0.5;
            if (survival(middle, eigenvalues).pValue() > survivalProbability)
                low = middle;
            else
                high = middle;
        }
        return high;
    }

    private static Double imhof(double statistic, double[] eigenvalues) {
        double maximum = Arrays.stream(eigenvalues).max().orElse(0);
        double blockWidth = 1.0 / maximum;
        double integral = 0;
        int quiet = 0;
        for (int block = 0; block < MAXIMUM_BLOCKS; block++) {
            double start = block * blockWidth;
            double contribution = simpson(start, start + blockWidth,
                statistic, eigenvalues);
            if (!Double.isFinite(contribution)) return null;
            integral += contribution;
            quiet = Math.abs(contribution) < ABSOLUTE_TOLERANCE
                ? quiet + 1 : 0;
            if (quiet >= REQUIRED_QUIET_BLOCKS) {
                double result = 0.5 + integral / Math.PI;
                return result >= -1e-7 && result <= 1 + 1e-7
                    ? Math.max(0, Math.min(1, result)) : null;
            }
        }
        return null;
    }

    private static double simpson(
            double start, double end, double statistic,
            double[] eigenvalues) {
        double step = (end - start) / INTERVALS_PER_BLOCK;
        double total = integrand(start, statistic, eigenvalues)
            + integrand(end, statistic, eigenvalues);
        for (int index = 1; index < INTERVALS_PER_BLOCK; index++)
            total += (index % 2 == 0 ? 2 : 4)
                * integrand(start + index * step, statistic, eigenvalues);
        return total * step / 3;
    }

    private static double integrand(
            double value, double statistic, double[] eigenvalues) {
        if (value == 0) return Arrays.stream(eigenvalues).sum() - statistic;
        double angle = -statistic * value;
        double logDenominator = 0;
        for (double eigenvalue : eigenvalues) {
            double term = 2 * eigenvalue * value;
            angle += 0.5 * Math.atan(term);
            logDenominator += 0.25 * Math.log1p(term * term);
        }
        return Math.sin(angle) * Math.exp(-logDenominator) / value;
    }

    private static double probability(double value) {
        return Math.max(Double.MIN_VALUE, Math.min(1, value));
    }
}
