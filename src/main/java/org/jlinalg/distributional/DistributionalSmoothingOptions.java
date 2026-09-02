/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

/** Outer AIC smoothing controls for multi-predictor distributional models. */
public record DistributionalSmoothingOptions(
        int maximumSweeps,
        double initialLogStep,
        double minimumLogSmoothing,
        double maximumLogSmoothing,
        double tolerance) {
    public DistributionalSmoothingOptions {
        if (maximumSweeps < 1 || !(initialLogStep > 0.0)
                || !Double.isFinite(initialLogStep)
                || !Double.isFinite(minimumLogSmoothing)
                || !(maximumLogSmoothing > minimumLogSmoothing)
                || !(tolerance > 0.0) || !Double.isFinite(tolerance)) {
            throw new IllegalArgumentException("invalid distributional smoothing controls");
        }
    }
    public static DistributionalSmoothingOptions defaults() {
        return new DistributionalSmoothingOptions(15, 2.0, -14.0, 18.0, 1e-3);
    }
}
