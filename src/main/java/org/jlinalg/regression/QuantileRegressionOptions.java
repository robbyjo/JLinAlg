/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

/** Controls for deterministic smoothed pinball-loss quantile regression. */
public record QuantileRegressionOptions(int maximumIterations, double relativeTolerance,
                                        double initialStep, double smoothing) {
    public QuantileRegressionOptions {
        if (maximumIterations < 1 || !(relativeTolerance > 0)
                || !(initialStep > 0) || !(smoothing > 0)
                || !Double.isFinite(relativeTolerance) || !Double.isFinite(initialStep)
                || !Double.isFinite(smoothing))
            throw new IllegalArgumentException("quantile controls are invalid");
    }

    public static QuantileRegressionOptions defaults() {
        return new QuantileRegressionOptions(2_000, 1e-8, 0.05, 1e-3);
    }
}
