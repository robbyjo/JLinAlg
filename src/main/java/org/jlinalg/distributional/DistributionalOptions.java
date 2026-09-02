/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

/** Immutable controls for block Fisher-scoring distributional fits. */
public record DistributionalOptions(
        int maximumIterations,
        double relativeTolerance,
        double maximumStep,
        double confidenceLevel) {
    public DistributionalOptions {
        if (maximumIterations < 1) {
            throw new IllegalArgumentException("maximumIterations must be positive");
        }
        if (!(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)) {
            throw new IllegalArgumentException(
                "relativeTolerance must be finite and positive");
        }
        if (!(maximumStep > 0.0) || !Double.isFinite(maximumStep)) {
            throw new IllegalArgumentException(
                "maximumStep must be finite and positive");
        }
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException(
                "confidenceLevel must lie strictly between zero and one");
        }
    }

    public static DistributionalOptions defaults() {
        return new DistributionalOptions(100, 1e-8, 3.0, 0.95);
    }
}
