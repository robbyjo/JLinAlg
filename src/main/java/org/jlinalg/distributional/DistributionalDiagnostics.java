/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import jdistlib.Beta;
import jdistlib.Gamma;
import jdistlib.NegBinomial;
import jdistlib.Normal;
import jdistlib.Poisson;

/** Quantile residuals and centiles for built-in distributional families. */
public final class DistributionalDiagnostics {
    private static final double EPSILON = 1e-12;

    private DistributionalDiagnostics() { }

    /**
     * Returns normal quantile residuals. Discrete families use the deterministic
     * midpoint of the probability interval; callers needing randomized Dunn-
     * Smyth residuals can randomize uniformly between the same CDF bounds.
     */
    public static double[] quantileResiduals(
            double[] response,
            DistributionalResult fit) {
        if (response == null || fit == null
                || response.length != fit.parameters().get(0).fittedValues().length) {
            throw new IllegalArgumentException("response must match fitted observations");
        }
        double[] result = new double[response.length];
        for (int row = 0; row < response.length; row++) {
            Bounds bounds = probabilityBounds(response[row], fit, row);
            double probability = 0.5 * (bounds.lower() + bounds.upper());
            probability = Math.max(EPSILON, Math.min(1.0 - EPSILON, probability));
            result[row] = Normal.quantile(probability, 0.0, 1.0, true, false);
        }
        return result;
    }

    /** Gaussian centile prediction from fitted location and scale. */
    public static double[] gaussianCentile(
            DistributionalResult fit, double probability) {
        if (!"gaussian-location-scale".equals(fit.family())
                || !(probability > 0.0 && probability < 1.0)) {
            throw new IllegalArgumentException(
                "Gaussian location-scale fit and probability in (0,1) are required");
        }
        double[] mean = fit.parameter("mu").fittedValues();
        double[] scale = fit.parameter("sigma").fittedValues();
        double quantile = Normal.quantile(probability, 0.0, 1.0, true, false);
        double[] result = new double[mean.length];
        for (int row = 0; row < mean.length; row++) {
            result[row] = mean[row] + scale[row] * quantile;
        }
        return result;
    }

    private static Bounds probabilityBounds(
            double response, DistributionalResult fit, int row) {
        return switch (fit.family()) {
            case "gaussian-location-scale" -> {
                double mean = value(fit, "mu", row);
                double scale = value(fit, "sigma", row);
                double probability = Normal.cumulative(
                    response, mean, scale, true, false);
                yield new Bounds(probability, probability);
            }
            case "gamma-mean-shape" -> {
                double mean = value(fit, "mu", row);
                double shape = value(fit, "shape", row);
                double probability = Gamma.cumulative(
                    response, shape, mean / shape, true, false);
                yield new Bounds(probability, probability);
            }
            case "beta-mean-precision" -> {
                double mean = value(fit, "mu", row);
                double precision = value(fit, "precision", row);
                double probability = Beta.cumulative(response,
                    mean * precision, (1.0 - mean) * precision, true, false);
                yield new Bounds(probability, probability);
            }
            case "negative-binomial-mean-size" -> {
                double mean = value(fit, "mu", row);
                double size = value(fit, "size", row);
                yield new Bounds(
                    response == 0.0 ? 0.0 : NegBinomial.cumulative_mu(
                        response - 1.0, size, mean, true, false),
                    NegBinomial.cumulative_mu(response, size, mean, true, false));
            }
            case "zero-inflated-poisson" -> {
                double mean = value(fit, "mu", row);
                double zero = value(fit, "zeroProbability", row);
                double lower = response == 0.0 ? 0.0 : zero
                    + (1.0 - zero) * Poisson.cumulative(
                        response - 1.0, mean, true, false);
                double upper = zero + (1.0 - zero) * Poisson.cumulative(
                    response, mean, true, false);
                yield new Bounds(lower, upper);
            }
            case "hurdle-poisson" -> hurdleBounds(response, fit, row);
            default -> throw new IllegalArgumentException(
                "quantile residuals are not defined for family " + fit.family());
        };
    }

    private static Bounds hurdleBounds(
            double response, DistributionalResult fit, int row) {
        double mean = value(fit, "mu", row);
        double zero = value(fit, "zeroProbability", row);
        if (response == 0.0) return new Bounds(0.0, zero);
        double poissonZero = Math.exp(-mean);
        double denominator = 1.0 - poissonZero;
        double lowerPositive = (Poisson.cumulative(
            response - 1.0, mean, true, false) - poissonZero) / denominator;
        double upperPositive = (Poisson.cumulative(
            response, mean, true, false) - poissonZero) / denominator;
        return new Bounds(zero + (1.0 - zero) * lowerPositive,
            zero + (1.0 - zero) * upperPositive);
    }

    private static double value(
            DistributionalResult fit, String parameter, int row) {
        return fit.parameter(parameter).fittedValues()[row];
    }

    private record Bounds(double lower, double upper) { }
}
