/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

/** Controls for frequentist zero-inflated sparse Laplace mixed models. */
public record ZeroInflatedMixedOptions(
        int maximumOuterEvaluations,
        int maximumModeIterations,
        double relativeTolerance,
        double initialTrustRadius,
        double minimumVariance,
        double maximumVariance,
        double minimumSize,
        double maximumSize,
        double maximumAbsoluteCoefficient,
        double[] initialVariances) {
    public ZeroInflatedMixedOptions {
        if (maximumOuterEvaluations < 20 || maximumModeIterations < 1
                || !(relativeTolerance > 0.0)
                || !Double.isFinite(relativeTolerance)
                || !(initialTrustRadius > 0.0)
                || !Double.isFinite(initialTrustRadius)
                || !(minimumVariance > 0.0)
                || !(maximumVariance > minimumVariance)
                || !(minimumSize > 0.0)
                || !(maximumSize > minimumSize)
                || !(maximumAbsoluteCoefficient > 0.0)
                || !Double.isFinite(maximumAbsoluteCoefficient)) {
            throw new IllegalArgumentException(
                "invalid zero-inflated mixed-model controls");
        }
        initialVariances = initialVariances == null
            ? null : initialVariances.clone();
        if (initialVariances != null) {
            for (double value : initialVariances) {
                if (!(value > 0.0) || !Double.isFinite(value)) {
                    throw new IllegalArgumentException(
                        "initial variances must be finite and positive");
                }
            }
        }
    }

    @Override public double[] initialVariances() {
        return initialVariances == null ? null : initialVariances.clone();
    }

    public static ZeroInflatedMixedOptions defaults() {
        return new ZeroInflatedMixedOptions(
            500, 80, 1e-7, 0.5,
            1e-8, 1e8, 1e-5, 1e7, 30.0, null);
    }
}
