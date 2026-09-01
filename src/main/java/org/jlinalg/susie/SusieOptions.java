/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

/** Controls for iterative Bayesian stepwise selection. */
public record SusieOptions(
        int effects,
        int maximumIterations,
        double convergenceTolerance,
        double priorVariance,
        boolean estimateResidualVariance,
        double credibleSetCoverage,
        double minimumCredibleSetPurity) {
    public SusieOptions {
        if (effects < 1 || maximumIterations < 1
                || !(convergenceTolerance > 0.0)
                || !(priorVariance > 0.0)
                || !(credibleSetCoverage > 0.0 && credibleSetCoverage < 1.0)
                || !(minimumCredibleSetPurity >= 0.0 && minimumCredibleSetPurity <= 1.0)) {
            throw new IllegalArgumentException("invalid SuSiE options");
        }
    }
    public static SusieOptions defaults() {
        return new SusieOptions(10, 200, 1e-6, 0.2, true, 0.95, 0.5);
    }
}
