/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

/** Controls for first-order Laplace marginal-likelihood GLMM fitting. */
public record GlmmLaplaceOptions(
        int maximumOuterIterations,
        int maximumModeIterations,
        double relativeTolerance,
        double initialLogVarianceStep,
        double minimumVariance,
        double maximumVariance,
        double[] initialVariances) {
    public GlmmLaplaceOptions {
        if (maximumOuterIterations < 1 || maximumModeIterations < 1
                || !(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)
                || !(initialLogVarianceStep > 0.0)
                || !Double.isFinite(initialLogVarianceStep)
                || !(minimumVariance > 0.0)
                || !(maximumVariance > minimumVariance)
                || !Double.isFinite(maximumVariance)) {
            throw new IllegalArgumentException("invalid Laplace GLMM controls");
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
    public static GlmmLaplaceOptions defaults() {
        return new GlmmLaplaceOptions(
            30, 80, 1e-7, 1.5, 1e-8, 1e8, null);
    }
}
