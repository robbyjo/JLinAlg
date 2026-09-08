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
        return jdistlib.NegBinomial.density_mu(response, size, mean, true);
    }
    @Override public void derivatives(
            double response, double[] parameters,
            double[] score, double[] information) {
        double mean = parameters[0];
        double size = parameters[1];
        double total = mean + size;
        score[0] = size * (response - mean) / total;
        score[1] = sizeScore(response, mean, size);
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
    private static double sizeScore(double y, double mu, double size) {
        if (Math.max(y, mu) / size < 1e-4) {
            double first = .5 * (y - (y - mu) * (y - mu));
            double second = y*(y-1)*(2*y-1)/6 + 2*mu*mu*mu/3 - y*mu*mu;
            return first / size + second / size / size;
        }
        return size * (SpecialFunctions.digamma(y + size) - SpecialFunctions.digamma(size)
            - Math.log1p(mu / size) + (mu - y) / (size + mu));
    }

    /** Actual observed information, distinct from the positive scoring metric. */
    void observedInformation(double y, double[] parameters, double[] information) {
        double mu = parameters[0], size = parameters[1], total = mu + size;
        information[0] = (size / total) * (mu / total) * (size + y);
        information[1] = information[2] = (size / total) * (mu / total) * (mu - y);
        if (Math.max(y, mu) / size < 1e-4) {
            double first = .5 * (y - (y - mu) * (y - mu));
            double second = y*(y-1)*(2*y-1)/6 + 2*mu*mu*mu/3 - y*mu*mu;
            information[3] = first / size + 2 * second / size / size;
        } else {
            information[3] = size * size * (SpecialFunctions.trigamma(size)
                - SpecialFunctions.trigamma(y + size) - mu / size / total
                + (mu - y) / total / total) - sizeScore(y, mu, size);
        }
    }
    private static void check(int parameter) {
        if (parameter < 0 || parameter >= 2) {
            throw new IllegalArgumentException("unknown NB parameter");
        }
    }
}
