/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Controls for deterministic softmax multinomial regression. */
public record MultinomialOptions(int maximumIterations, double relativeTolerance,
                                 double initialStep) {
    public MultinomialOptions {
        if (maximumIterations < 1 || !(relativeTolerance > 0)
                || !(initialStep > 0) || !Double.isFinite(relativeTolerance)
                || !Double.isFinite(initialStep))
            throw new IllegalArgumentException("multinomial controls are invalid");
    }
    public static MultinomialOptions defaults() { return new MultinomialOptions(500, 1e-8, 1.0); }
}
