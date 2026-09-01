/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;
import java.util.Arrays;

/** Two-component valid/invalid instrument likelihood over a causal-effect grid. */
public final class ContaminationMixture {
    private ContaminationMixture() { }

    public static ContaminationMixtureResult fit(
            List<HarmonizedInstrument> instruments, int gridPoints) {
        List<HarmonizedInstrument> values = MendelianRandomization.validated(instruments, 3);
        if (gridPoints < 101) throw new IllegalArgumentException("gridPoints must be at least 101");
        double minimum = Double.POSITIVE_INFINITY;
        double maximum = Double.NEGATIVE_INFINITY;
        for (HarmonizedInstrument value : values) {
            double ratio = value.outcomeEffect() / value.exposureEffect();
            minimum = Math.min(minimum, ratio);
            maximum = Math.max(maximum, ratio);
        }
        double span = Math.max(0.1, maximum - minimum);
        minimum -= span;
        maximum += span;
        double bestBeta = 0.0;
        double bestPi = 0.5;
        double best = Double.NEGATIVE_INFINITY;
        double[] profile = new double[gridPoints];
        Arrays.fill(profile, Double.NEGATIVE_INFINITY);
        for (int point = 0; point < gridPoints; point++) {
            double beta = minimum + (maximum - minimum) * point / (gridPoints - 1.0);
            for (int piIndex = 1; piIndex < 20; piIndex++) {
                double pi = piIndex / 20.0;
                double likelihood = likelihood(values, beta, pi);
                if (likelihood > profile[point]) profile[point] = likelihood;
                if (likelihood > best) {
                    best = likelihood;
                    bestBeta = beta;
                    bestPi = pi;
                }
            }
        }
        double step = (maximum - minimum) / (gridPoints - 1.0);
        int bestIndex = (int) Math.round((bestBeta - minimum) / step);
        double curvature = bestIndex > 0 && bestIndex < gridPoints - 1
            ? -(profile[bestIndex + 1] - 2.0 * profile[bestIndex]
                + profile[bestIndex - 1]) / (step * step) : Double.NaN;
        double se = curvature > 0.0 ? 1.0 / Math.sqrt(curvature) : step;
        MrEstimate estimate = MendelianRandomization.estimate(
            MrMethod.CONTAMINATION_MIXTURE, bestBeta, se, 0.95,
            Double.NaN, 0, Double.NaN, values.size());
        return new ContaminationMixtureResult(estimate, bestPi, best, gridPoints);
    }

    private static double likelihood(
            List<HarmonizedInstrument> values, double beta, double pi) {
        double result = 0.0;
        for (HarmonizedInstrument value : values) {
            double residual = value.outcomeEffect() - beta * value.exposureEffect();
            double variance = value.outcomeStandardError() * value.outcomeStandardError()
                + beta * beta * value.exposureStandardError() * value.exposureStandardError();
            double valid = normalDensity(residual, variance);
            double invalid = normalDensity(residual, variance + 0.01);
            result += Math.log(Math.max(1e-300, pi * valid + (1.0 - pi) * invalid));
        }
        return result;
    }

    private static double normalDensity(double value, double variance) {
        return Math.exp(-0.5 * value * value / variance)
            / Math.sqrt(2.0 * Math.PI * variance);
    }
}
