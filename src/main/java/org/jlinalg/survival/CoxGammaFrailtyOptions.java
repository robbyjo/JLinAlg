/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Controls for Laplace-profiled shared gamma frailty. */
public record CoxGammaFrailtyOptions(
        CoxOptions coxOptions,
        double initialVariance,
        double minimumVariance,
        double maximumVariance,
        int maximumVarianceIterations,
        double logVarianceTolerance) {
    public CoxGammaFrailtyOptions {
        if (coxOptions == null || !(minimumVariance > 0.0)
                || !(maximumVariance > minimumVariance)
                || initialVariance < minimumVariance
                || initialVariance > maximumVariance
                || maximumVarianceIterations < 1
                || !(logVarianceTolerance > 0.0))
            throw new IllegalArgumentException(
                "invalid gamma-frailty controls");
    }

    public static CoxGammaFrailtyOptions defaults() {
        return new CoxGammaFrailtyOptions(CoxOptions.defaults(), 0.1,
            1e-6, 100.0, 30, 1e-4);
    }
}
