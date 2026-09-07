/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

/** Distribution, link, variance, and deviance contract for IRLS. */
public interface GlmFamily {
    String name();
    void validateResponse(double response, double priorWeight);
    double initialMean(double response);
    double link(double mean);
    double inverseLink(double linearPredictor);
    double meanDerivative(double linearPredictor);
    double variance(double mean);
    double unitDeviance(double response, double mean);
    double logLikelihood(
        double response, double mean, double priorWeight, double dispersion);
    boolean fixedDispersion();

    /** Returns the Fisher information for the linear predictor. */
    default double workingWeight(
            double response, double linearPredictor, double mean,
            double priorWeight) {
        double derivative = meanDerivative(linearPredictor);
        return priorWeight * derivative * derivative / variance(mean);
    }

    /** Returns the offset-free Fisher-scoring working response. */
    default double workingResponse(
            double response, double linearPredictor, double mean,
            double priorWeight, double offset) {
        return linearPredictor - offset
            + (response - mean) / meanDerivative(linearPredictor);
    }
}
