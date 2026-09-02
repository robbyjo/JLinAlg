/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import java.util.List;

/** Beta distribution parameterized by modeled mean and positive precision. */
public final class BetaMeanPrecisionFamily implements DistributionalFamily {
    private static final double EPSILON = 1e-10;

    @Override public String name() { return "beta-mean-precision"; }
    @Override public int parameterCount() { return 2; }
    @Override public List<String> parameterNames() {
        return List.of("mu", "precision");
    }
    @Override public void validateResponse(double response) {
        if (!(response > 0.0 && response < 1.0)
                || !Double.isFinite(response)) {
            throw new IllegalArgumentException(
                "Beta responses must lie strictly between zero and one");
        }
    }
    @Override public double[] initialParameters(double[] response) {
        double mean = Arrays.stream(response).average().orElseThrow();
        double variance = 0.0;
        for (double value : response) {
            double difference = value - mean;
            variance += difference * difference;
        }
        variance /= Math.max(1, response.length - 1);
        double precision = variance == 0.0 ? 100.0
            : mean * (1.0 - mean) / variance - 1.0;
        return new double[] {
            clamp(mean, EPSILON, 1.0 - EPSILON),
            Math.max(EPSILON, precision)};
    }
    @Override public double link(int parameter, double value) {
        check(parameter);
        if (parameter == 0) {
            if (!(value > 0.0 && value < 1.0)) {
                throw new IllegalArgumentException("Beta mean must lie in (0,1)");
            }
            return Math.log(value / (1.0 - value));
        }
        if (!(value > 0.0)) {
            throw new IllegalArgumentException("Beta precision must be positive");
        }
        return Math.log(value);
    }
    @Override public double inverseLink(int parameter, double predictor) {
        check(parameter);
        if (parameter == 0) {
            double probability = predictor >= 0.0
                ? 1.0 / (1.0 + Math.exp(-predictor))
                : Math.exp(predictor) / (1.0 + Math.exp(predictor));
            return clamp(probability, EPSILON, 1.0 - EPSILON);
        }
        return Math.max(EPSILON, Math.exp(Math.min(350.0, predictor)));
    }
    @Override public double logLikelihood(
            double response, double[] parameters) {
        double mean = parameters[0];
        double precision = parameters[1];
        double alpha = mean * precision;
        double beta = (1.0 - mean) * precision;
        return SpecialFunctions.logGamma(precision)
            - SpecialFunctions.logGamma(alpha)
            - SpecialFunctions.logGamma(beta)
            + (alpha - 1.0) * Math.log(response)
            + (beta - 1.0) * Math.log1p(-response);
    }
    @Override public void derivatives(
            double response,
            double[] parameters,
            double[] score,
            double[] information) {
        double mean = parameters[0];
        double precision = parameters[1];
        double alpha = mean * precision;
        double beta = (1.0 - mean) * precision;
        double meanDerivative = mean * (1.0 - mean);
        double logitResponse = Math.log(response) - Math.log1p(-response);
        double meanScore = precision * (-SpecialFunctions.digamma(alpha)
            + SpecialFunctions.digamma(beta) + logitResponse);
        double precisionScore = SpecialFunctions.digamma(precision)
            - mean * SpecialFunctions.digamma(alpha)
            - (1.0 - mean) * SpecialFunctions.digamma(beta)
            + mean * Math.log(response)
            + (1.0 - mean) * Math.log1p(-response);
        score[0] = meanDerivative * meanScore;
        score[1] = precision * precisionScore;

        double trigammaAlpha = SpecialFunctions.trigamma(alpha);
        double trigammaBeta = SpecialFunctions.trigamma(beta);
        double informationMean = precision * precision
            * (trigammaAlpha + trigammaBeta);
        double informationCross = precision
            * (mean * trigammaAlpha
                - (1.0 - mean) * trigammaBeta);
        double informationPrecision = -SpecialFunctions.trigamma(precision)
            + mean * mean * trigammaAlpha
            + (1.0 - mean) * (1.0 - mean) * trigammaBeta;
        information[0] = meanDerivative * meanDerivative * informationMean;
        information[1] = meanDerivative * precision * informationCross;
        information[2] = information[1];
        information[3] = precision * precision * informationPrecision;
    }
    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }
    private static void check(int parameter) {
        if (parameter < 0 || parameter >= 2) {
            throw new IllegalArgumentException("unknown Beta parameter");
        }
    }
}
