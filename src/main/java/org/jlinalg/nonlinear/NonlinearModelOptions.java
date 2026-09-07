/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

/** Deterministic Gauss-Newton controls for nonlinear fixed and mixed models. */
public record NonlinearModelOptions(int maximumIterations,
        double parameterTolerance, double objectiveTolerance,
        double initialStep, int maximumStepHalvings) {
    public NonlinearModelOptions {
        if (maximumIterations < 1 || !(parameterTolerance > 0)
                || !(objectiveTolerance > 0) || !(initialStep > 0)
                || maximumStepHalvings < 1
                || !Double.isFinite(parameterTolerance)
                || !Double.isFinite(objectiveTolerance)
                || !Double.isFinite(initialStep))
            throw new IllegalArgumentException("nonlinear optimization controls are invalid");
    }

    public static NonlinearModelOptions defaults() {
        return new NonlinearModelOptions(100, 1e-8, 1e-10, 1.0, 20);
    }
}
