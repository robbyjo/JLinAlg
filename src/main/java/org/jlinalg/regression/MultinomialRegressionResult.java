/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Multinomial logistic coefficients and fitted class probabilities. */
public record MultinomialRegressionResult(
        double[] coefficients, double[] probabilities, double logLikelihood,
        int observations, int predictorCount, int classCount,
        int iterations, boolean converged) {
    public MultinomialRegressionResult { coefficients = coefficients.clone(); probabilities = probabilities.clone(); }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] probabilities() { return probabilities.clone(); }
}
