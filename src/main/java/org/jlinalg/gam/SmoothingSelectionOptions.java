/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

/** Controls for multi-penalty Gaussian smoothing selection. */
public record SmoothingSelectionOptions(
        SmoothingCriterion criterion,
        int maximumSweeps,
        double initialLogStep,
        double minimumLogSmoothing,
        double maximumLogSmoothing,
        double tolerance,
        double knownScale) {
    public SmoothingSelectionOptions {
        if (criterion == null || maximumSweeps < 1
                || !(initialLogStep > 0.0) || !Double.isFinite(initialLogStep)
                || !Double.isFinite(minimumLogSmoothing)
                || !Double.isFinite(maximumLogSmoothing)
                || !(maximumLogSmoothing > minimumLogSmoothing)
                || !(tolerance > 0.0) || !Double.isFinite(tolerance)
                || (criterion == SmoothingCriterion.UBRE
                    && (!(knownScale > 0.0) || !Double.isFinite(knownScale)))) {
            throw new IllegalArgumentException("invalid smoothing-selection controls");
        }
    }

    /** Default GCV controls spanning smoothing parameters from exp(-16) to exp(20). */
    public static SmoothingSelectionOptions gcv() {
        return new SmoothingSelectionOptions(
            SmoothingCriterion.GCV, 20, 3.0, -16.0, 20.0, 1e-4, Double.NaN);
    }

    /** UBRE controls for a known residual variance. */
    public static SmoothingSelectionOptions ubre(double scale) {
        return new SmoothingSelectionOptions(
            SmoothingCriterion.UBRE, 20, 3.0, -16.0, 20.0, 1e-4, scale);
    }
}
