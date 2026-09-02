/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import java.util.List;

/** NB2 counts with separate additive predictors for mean and size. */
public final class NegativeBinomialMeanDispersionFamily
        implements DistributionalFamily {
    private static final double MINIMUM = 1e-10;

    @Override public String name() { return "negative-binomial-mean-size"; }
    @Override public int parameterCount() { return 2; }
    @Override public List<String> parameterNames() { return List.of("mu", "size"); }
    @Override public void validateResponse(double response) {
        if (!Double.isFinite(response) || response < 0.0
                || response != Math.rint(response)) {
            throw new IllegalArgumentException(
                "negative-binomial responses must be nonnegative integers");
        }
    }
    @Override public double[] initialParameters(double[] response) {
        double mean = Arrays.stream(response).average().orElseThrow();
        double variance = 0.0;
        for (double value : response) variance += (value - mean) * (value - mean);
        variance /= Math.max(1, response.length - 1);
        double size = variance > mean
            ? mean * mean / Math.max(MINIMUM, variance - mean) : 100.0;
        return new double[] {Math.max(MINIMUM, mean), Math.max(MINIMUM, size)};
    }
    @Override public double link(int parameter, double value) {
        check(parameter);
        if (!(value > 0.0)) {
            throw new IllegalArgumentException("NB parameters must be positive");
        }
        return Math.log(value);
    }
    @Override public double inverseLink(int parameter, double predictor) {
        check(parameter);
        return Math.max(MINIMUM, Math.exp(Math.min(350.0, predictor)));
    }
    @Override public double logLikelihood(double response, double[] parameters) {
        double mean = parameters[0];
        double size = parameters[1];
        return SpecialFunctions.logGamma(response + size)
            - SpecialFunctions.logGamma(size)
            - SpecialFunctions.logGamma(response + 1.0)
            + size * (Math.log(size) - Math.log(size + mean))
            + response * (Math.log(mean) - Math.log(size + mean));
    }
    @Override public void derivatives(
            double response, double[] parameters,
            double[] score, double[] information) {
        double mean = parameters[0];
        double size = parameters[1];
        double total = mean + size;
        score[0] = size * (response - mean) / total;
        score[1] = size * (SpecialFunctions.digamma(response + size)
            - SpecialFunctions.digamma(size) + Math.log(size) + 1.0
            - Math.log(total) - (response + size) / total);
        // Mean and size are orthogonal in expectation; an OPG contribution
        // stabilizes the size update when the expected trigamma term is small.
        information[0] = Math.max(MINIMUM, size * mean / total);
        information[1] = 0.0;
        information[2] = 0.0;
        double curvature = size * size * Math.max(MINIMUM,
            SpecialFunctions.trigamma(size)
                - SpecialFunctions.trigamma(response + size)
                - 1.0 / size + 1.0 / total
                + (mean - response) / (total * total));
        information[3] = Math.max(1e-6, curvature + 0.05 * score[1] * score[1]);
    }
    private static void check(int parameter) {
        if (parameter < 0 || parameter >= 2) {
            throw new IllegalArgumentException("unknown NB parameter");
        }
    }
}
