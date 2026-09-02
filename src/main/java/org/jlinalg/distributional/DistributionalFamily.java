/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.List;

/** Likelihood, links, scores, and Fisher information for vector predictors. */
public interface DistributionalFamily {
    String name();
    int parameterCount();
    List<String> parameterNames();
    void validateResponse(double response);
    double[] initialParameters(double[] response);
    double link(int parameter, double value);
    double inverseLink(int parameter, double predictor);
    double logLikelihood(double response, double[] parameters);

    /**
     * Converts linked vector predictors to category probabilities when this is
     * a discrete vector-response family.
     */
    default double[] categoryProbabilities(double[] parameters) {
        throw new UnsupportedOperationException(
            name() + " does not define category probabilities");
    }

    /**
     * Fills the score with derivatives with respect to linked predictors and
     * fills the row-major expected negative Hessian (Fisher information).
     */
    void derivatives(
        double response,
        double[] parameters,
        double[] score,
        double[] information);
}
