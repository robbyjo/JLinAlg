/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import jdistlib.ChiSquare;

/** Likelihood-based comparison of nested distributional models. */
public record DistributionalModelComparison(
        double likelihoodRatio,
        int degreesOfFreedom,
        double pValue,
        double smallerAic,
        double largerAic) {
    /** Compares a nested smaller model to a larger model using total coefficient counts. */
    public static DistributionalModelComparison compare(
            DistributionalResult smaller,
            DistributionalResult larger) {
        if (smaller == null || larger == null
                || !smaller.family().equals(larger.family())) {
            throw new IllegalArgumentException(
                "nested fits from the same family are required");
        }
        int smallerParameters = coefficientCount(smaller);
        int largerParameters = coefficientCount(larger);
        int degrees = largerParameters - smallerParameters;
        if (degrees < 1) {
            throw new IllegalArgumentException(
                "larger model must contain more coefficients");
        }
        double statistic = Math.max(0.0,
            2.0 * (larger.logLikelihood() - smaller.logLikelihood()));
        double p = ChiSquare.cumulative(statistic, degrees, false, false);
        return new DistributionalModelComparison(statistic, degrees, p,
            -2.0 * smaller.logLikelihood() + 2.0 * smallerParameters,
            -2.0 * larger.logLikelihood() + 2.0 * largerParameters);
    }
    private static int coefficientCount(DistributionalResult fit) {
        int result = 0;
        for (DistributionalParameterResult parameter : fit.parameters()) {
            result += parameter.coefficients().length;
        }
        return result;
    }
}
