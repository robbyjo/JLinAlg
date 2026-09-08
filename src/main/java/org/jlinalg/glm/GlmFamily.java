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

    /** Predictor-aware variance; families may avoid rounded inverse-link tails. */
    default double varianceAtPredictor(double predictor, double mean) { return variance(mean); }
    /** Predictor-aware response residual, including small complement probabilities. */
    default double residualAtPredictor(double response, double predictor, double mean) {
        return response - mean;
    }
    /** Predictor-aware deviance, without forcing tail probabilities through rounded means. */
    default double unitDevianceAtPredictor(double response, double predictor, double mean) {
        return unitDeviance(response, mean);
    }
    /** Predictor-aware density evaluation; default preserves custom-family contracts. */
    default double logLikelihoodAtPredictor(double response, double predictor, double mean,
            double priorWeight, double dispersion) {
        return logLikelihood(response, mean, priorWeight, dispersion);
    }

    /** Returns the Fisher information for the linear predictor. */
    default double workingWeight(
            double response, double linearPredictor, double mean,
            double priorWeight) {
        double derivative = meanDerivative(linearPredictor);
        double standardized = derivative / Math.sqrt(varianceAtPredictor(linearPredictor, mean));
        return priorWeight * standardized * standardized;
    }

    /** Returns the offset-free Fisher-scoring working response. */
    default double workingResponse(
            double response, double linearPredictor, double mean,
            double priorWeight, double offset) {
        return linearPredictor - offset
            + residualAtPredictor(response, linearPredictor, mean) / meanDerivative(linearPredictor);
    }
}
