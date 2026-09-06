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
        double[] initialVariances,
        ZeroInflatedOuterOptimizer outerOptimizer,
        int maximumGradientThreads) {
    /** Backward-compatible controls using automatic optimizer selection. */
    public ZeroInflatedMixedOptions(
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
        this(maximumOuterEvaluations, maximumModeIterations,
            relativeTolerance, initialTrustRadius, minimumVariance,
            maximumVariance, minimumSize, maximumSize,
            maximumAbsoluteCoefficient, initialVariances,
            ZeroInflatedOuterOptimizer.AUTO, 4);
    }

    /** Controls with explicit outer optimizer and up to four gradient threads. */
    public ZeroInflatedMixedOptions(
            int maximumOuterEvaluations,
            int maximumModeIterations,
            double relativeTolerance,
            double initialTrustRadius,
            double minimumVariance,
            double maximumVariance,
            double minimumSize,
            double maximumSize,
            double maximumAbsoluteCoefficient,
            double[] initialVariances,
            ZeroInflatedOuterOptimizer outerOptimizer) {
        this(maximumOuterEvaluations, maximumModeIterations,
            relativeTolerance, initialTrustRadius, minimumVariance,
            maximumVariance, minimumSize, maximumSize,
            maximumAbsoluteCoefficient, initialVariances, outerOptimizer, 4);
    }

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
                || !Double.isFinite(maximumAbsoluteCoefficient)
                || maximumGradientThreads < 1) {
            throw new IllegalArgumentException(
                "invalid zero-inflated mixed-model controls");
        }
        initialVariances = initialVariances == null
            ? null : initialVariances.clone();
        outerOptimizer = outerOptimizer == null
            ? ZeroInflatedOuterOptimizer.AUTO : outerOptimizer;
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
            1e-8, 1e8, 1e-5, 1e7, 30.0, null,
            ZeroInflatedOuterOptimizer.AUTO, 4);
    }
}
