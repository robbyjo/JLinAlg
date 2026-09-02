/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import java.util.List;

/** Poisson positive counts with a separate logit model for the zero hurdle. */
public final class HurdlePoissonFamily implements DistributionalFamily {
    private static final double EPSILON = 1e-10;

    @Override public String name() { return "hurdle-poisson"; }
    @Override public int parameterCount() { return 2; }
    @Override public List<String> parameterNames() {
        return List.of("mu", "zeroProbability");
    }
    @Override public void validateResponse(double response) {
        if (!Double.isFinite(response) || response < 0.0
                || response != Math.rint(response)) {
            throw new IllegalArgumentException(
                "hurdle-Poisson responses must be nonnegative integers");
        }
    }
    @Override public double[] initialParameters(double[] response) {
        double positiveSum = 0.0;
        int positives = 0;
        for (double value : response) {
            if (value > 0.0) {
                positives++;
                positiveSum += value;
            }
        }
        double zero = 1.0 - positives / (double) response.length;
        return new double[] {
            Math.max(0.1, positives == 0 ? 0.1 : positiveSum / positives),
            clamp(zero)};
    }
    @Override public double link(int parameter, double value) {
        check(parameter);
        if (parameter == 0) {
            if (!(value > 0.0)) throw new IllegalArgumentException("mu must be positive");
            return Math.log(value);
        }
        value = clamp(value);
        return Math.log(value / (1.0 - value));
    }
    @Override public double inverseLink(int parameter, double predictor) {
        check(parameter);
        if (parameter == 0) {
            return Math.max(EPSILON, Math.exp(Math.min(350.0, predictor)));
        }
        return clamp(predictor >= 0.0
            ? 1.0 / (1.0 + Math.exp(-predictor))
            : Math.exp(predictor) / (1.0 + Math.exp(predictor)));
    }
    @Override public double logLikelihood(double response, double[] parameters) {
        double mean = parameters[0];
        double zero = parameters[1];
        if (response == 0.0) return Math.log(zero);
        return Math.log1p(-zero) + response * Math.log(mean) - mean
            - SpecialFunctions.logGamma(response + 1.0)
            - Math.log(-Math.expm1(-mean));
    }
    @Override public void derivatives(
            double response, double[] parameters,
            double[] score, double[] information) {
        double mean = parameters[0];
        double zero = parameters[1];
        if (response == 0.0) {
            score[0] = 0.0;
            score[1] = 1.0 - zero;
        } else {
            score[0] = response - mean - mean / Math.expm1(mean);
            score[1] = -zero;
        }
        Arrays.fill(information, 0.0);
        information[3] = Math.max(EPSILON, zero * (1.0 - zero));
        double exponentialMinusOne = Math.expm1(mean);
        double derivative = (exponentialMinusOne
            - mean * Math.exp(mean))
            / (exponentialMinusOne * exponentialMinusOne);
        information[0] = Math.max(EPSILON,
            (1.0 - zero) * (mean + mean * derivative));
    }
    private static double clamp(double value) {
        return Math.max(EPSILON, Math.min(1.0 - EPSILON, value));
    }
    private static void check(int parameter) {
        if (parameter < 0 || parameter >= 2) {
            throw new IllegalArgumentException("unknown hurdle-Poisson parameter");
        }
    }
}
