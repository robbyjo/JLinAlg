/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.Arrays;
import java.util.List;

/** Gamma distribution parameterized by modeled mean and positive shape. */
public final class GammaMeanShapeFamily implements DistributionalFamily {
    private static final double MINIMUM = 1e-10;

    @Override public String name() { return "gamma-mean-shape"; }
    @Override public int parameterCount() { return 2; }
    @Override public List<String> parameterNames() {
        return List.of("mu", "shape");
    }
    @Override public void validateResponse(double response) {
        if (!(response > 0.0) || !Double.isFinite(response)) {
            throw new IllegalArgumentException(
                "Gamma responses must be finite and positive");
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
        double shape = variance == 0.0 ? 100.0 : mean * mean / variance;
        return new double[] {mean, Math.max(MINIMUM, shape)};
    }
    @Override public double link(int parameter, double value) {
        check(parameter);
        if (!(value > 0.0)) {
            throw new IllegalArgumentException(
                "Gamma mean and shape must be positive");
        }
        return Math.log(value);
    }
    @Override public double inverseLink(int parameter, double predictor) {
        check(parameter);
        return Math.max(MINIMUM, Math.exp(Math.min(350.0, predictor)));
    }
    @Override public double logLikelihood(
            double response, double[] parameters) {
        double mean = parameters[0];
        double shape = parameters[1];
        return (shape - 1.0) * Math.log(response)
            - shape * response / mean + shape * Math.log(shape)
            - shape * Math.log(mean) - SpecialFunctions.logGamma(shape);
    }
    @Override public void derivatives(
            double response,
            double[] parameters,
            double[] score,
            double[] information) {
        double mean = parameters[0];
        double shape = parameters[1];
        score[0] = shape * (response / mean - 1.0);
        double shapeScore = Math.log(response) - response / mean
            + Math.log(shape) + 1.0 - Math.log(mean)
            - SpecialFunctions.digamma(shape);
        score[1] = shape * shapeScore;
        Arrays.fill(information, 0.0);
        information[0] = shape;
        information[3] = Math.max(MINIMUM,
            shape * shape * SpecialFunctions.trigamma(shape) - shape);
    }
    private static void check(int parameter) {
        if (parameter < 0 || parameter >= 2) {
            throw new IllegalArgumentException("unknown Gamma parameter");
        }
    }
}
