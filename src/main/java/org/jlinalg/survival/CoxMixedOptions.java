/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Laplace-profile controls for Gaussian Cox frailty models. */
public record CoxMixedOptions(
        CoxOptions coxOptions,
        double[] initialVariances,
        int maximumVarianceIterations,
        double logVarianceTolerance,
        double minimumVariance,
        double maximumVariance) {
    public CoxMixedOptions {
        if (coxOptions == null || initialVariances == null
                || initialVariances.length == 0
                || maximumVarianceIterations < 1
                || !(logVarianceTolerance > 0)
                || !Double.isFinite(logVarianceTolerance)
                || !(minimumVariance > 0)
                || !(maximumVariance > minimumVariance)
                || !Double.isFinite(maximumVariance))
            throw new IllegalArgumentException("invalid mixed Cox options");
        initialVariances = initialVariances.clone();
        for (double value : initialVariances)
            if (!(value >= minimumVariance && value <= maximumVariance)
                    || !Double.isFinite(value))
                throw new IllegalArgumentException(
                    "initial frailty variances must lie within bounds");
    }
    @Override public double[] initialVariances() {
        return initialVariances.clone();
    }

    public static CoxMixedOptions defaults() {
        return new CoxMixedOptions(CoxOptions.defaults(),
            new double[] {0.5}, 30, 1e-4, 1e-8, 1e4);
    }
}
