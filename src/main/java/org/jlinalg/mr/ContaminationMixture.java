/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;
import java.util.Arrays;

/**
 * Two-component valid/invalid instrument likelihood over a causal-effect grid.
 * Invalid effects have an additional fixed outcome-effect variance of 0.01;
 * this scale-dependent model is not an implementation of R's MRCML or ConMix.
 */
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
        if (!Double.isFinite(minimum) || !Double.isFinite(maximum))
            throw new IllegalArgumentException("ratio grid exceeds numerical range");
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
        if (!(curvature > 0.0) || !Double.isFinite(curvature) || !Double.isFinite(best))
            throw new IllegalArgumentException("grid likelihood has no resolved interior curvature; increase resolution or revise the model");
        double se = 1.0 / Math.sqrt(curvature);
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
            double valid = Math.log(pi) + normalLogDensity(residual, variance);
            double invalid = Math.log1p(-pi) + normalLogDensity(residual, variance + 0.01);
            double maximum = Math.max(valid, invalid);
            if (maximum == Double.NEGATIVE_INFINITY) return maximum;
            result += maximum + Math.log1p(Math.exp(Math.min(valid, invalid) - maximum));
        }
        return result;
    }

    private static double normalLogDensity(double value, double variance) {
        double standardized = value / Math.sqrt(variance);
        return -0.5 * standardized * standardized
            - 0.5 * (Math.log(2.0 * Math.PI) + Math.log(variance));
    }
}
