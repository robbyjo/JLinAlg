/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

/** Immutable controls for maximum-likelihood beta regression. */
public record BetaRegressionOptions(
        int maximumIterations,
        double relativeTolerance,
        double maximumStep,
        BetaMeanLink meanLink,
        BetaPrecisionLink precisionLink) {
    public BetaRegressionOptions {
        if (maximumIterations < 1) {
            throw new IllegalArgumentException("maximumIterations must be positive");
        }
        if (!(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)) {
            throw new IllegalArgumentException(
                "relativeTolerance must be finite and positive");
        }
        if (!(maximumStep > 0.0) || !Double.isFinite(maximumStep)) {
            throw new IllegalArgumentException("maximumStep must be finite and positive");
        }
        if (meanLink == null || precisionLink == null) {
            throw new IllegalArgumentException("mean and precision links are required");
        }
    }

    /** Matches {@code betareg(y ~ x)}: logit mean and identity-linked constant precision. */
    public static BetaRegressionOptions constantPrecisionDefaults() {
        return new BetaRegressionOptions(
            100, 1e-8, 1_000.0, BetaMeanLink.LOGIT,
            BetaPrecisionLink.IDENTITY);
    }

    /** Matches {@code betareg(y ~ x | z)}: logit mean and log-linked precision. */
    public static BetaRegressionOptions variablePrecisionDefaults() {
        return new BetaRegressionOptions(
            200, 1e-8, 1_000.0, BetaMeanLink.LOGIT,
            BetaPrecisionLink.LOG);
    }
}
