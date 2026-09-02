/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import java.util.List;

/** Zero-inflated Poisson with modeled count mean and structural-zero probability. */
public final class ZeroInflatedPoissonFamily implements DistributionalFamily {
    private static final double EPSILON = 1e-10;

    @Override public String name() { return "zero-inflated-poisson"; }
    @Override public int parameterCount() { return 2; }
    @Override public List<String> parameterNames() {
        return List.of("mu", "zeroProbability");
    }
    @Override public void validateResponse(double response) {
        if (response < 0.0 || response != Math.rint(response)
                || !Double.isFinite(response)) {
            throw new IllegalArgumentException(
                "ZIP responses must be nonnegative integers");
        }
    }
    @Override public double[] initialParameters(double[] response) {
        double mean = Arrays.stream(response).average().orElseThrow();
        int zeros = 0;
        for (double value : response) if (value == 0.0) zeros++;
        double observedZero = zeros / (double) response.length;
        double poissonZero = Math.exp(-Math.max(EPSILON, mean));
        double inflation = Math.max(EPSILON,
            Math.min(1.0 - EPSILON,
                (observedZero - poissonZero) / (1.0 - poissonZero)));
        double countMean = Math.max(EPSILON, mean / (1.0 - inflation));
        return new double[] {countMean, inflation};
    }
    @Override public double link(int parameter, double value) {
        check(parameter);
        if (parameter == 0) {
            if (!(value > 0.0)) {
                throw new IllegalArgumentException("ZIP mean must be positive");
            }
            return Math.log(value);
        }
        if (!(value > 0.0 && value < 1.0)) {
            throw new IllegalArgumentException(
                "zero probability must lie in (0,1)");
        }
        return Math.log(value / (1.0 - value));
    }
    @Override public double inverseLink(int parameter, double predictor) {
        check(parameter);
        if (parameter == 0) {
            return Math.max(EPSILON, Math.exp(Math.min(350.0, predictor)));
        }
        double probability = predictor >= 0.0
            ? 1.0 / (1.0 + Math.exp(-predictor))
            : Math.exp(predictor) / (1.0 + Math.exp(predictor));
        return Math.max(EPSILON, Math.min(1.0 - EPSILON, probability));
    }
    @Override public double logLikelihood(
            double response, double[] parameters) {
        double mean = parameters[0];
        double inflation = parameters[1];
        if (response == 0.0) {
            return Math.log(inflation + (1.0 - inflation) * Math.exp(-mean));
        }
        return Math.log1p(-inflation) + response * Math.log(mean) - mean
            - SpecialFunctions.logGamma(response + 1.0);
    }
    @Override public void derivatives(
            double response,
            double[] parameters,
            double[] score,
            double[] information) {
        double mean = parameters[0];
        double inflation = parameters[1];
        double poissonZero = Math.exp(-mean);
        double zeroProbability = inflation + (1.0 - inflation) * poissonZero;
        double zeroCountPosterior =
            (1.0 - inflation) * poissonZero / zeroProbability;
        double zeroMeanScore = -zeroCountPosterior * mean;
        double zeroInflationScore = inflation * (1.0 - inflation)
            * (1.0 - poissonZero) / zeroProbability;
        if (response == 0.0) {
            score[0] = zeroMeanScore;
            score[1] = zeroInflationScore;
        } else {
            score[0] = response - mean;
            score[1] = -inflation;
        }

        double positiveMass = (1.0 - inflation) * (1.0 - poissonZero);
        information[0] = zeroProbability * zeroMeanScore * zeroMeanScore
            + (1.0 - inflation) * (mean - poissonZero * mean * mean);
        information[1] = zeroProbability
            * zeroMeanScore * zeroInflationScore
            - inflation * (1.0 - inflation) * mean * poissonZero;
        information[2] = information[1];
        information[3] = zeroProbability
            * zeroInflationScore * zeroInflationScore
            + positiveMass * inflation * inflation;
    }
    private static void check(int parameter) {
        if (parameter < 0 || parameter >= 2) {
            throw new IllegalArgumentException("unknown ZIP parameter");
        }
    }
}
