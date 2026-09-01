/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Newton and tie-handling controls for Cox models. */
public record CoxOptions(
        CoxTies ties,
        int maximumIterations,
        double relativeTolerance,
        double scoreTolerance,
        int maximumStepHalvings,
        double informationRidge,
        double confidenceLevel) {
    public CoxOptions {
        if (ties == null || maximumIterations < 1
                || !(relativeTolerance > 0) || !Double.isFinite(relativeTolerance)
                || !(scoreTolerance > 0) || !Double.isFinite(scoreTolerance)
                || maximumStepHalvings < 0
                || !(informationRidge >= 0) || !Double.isFinite(informationRidge)
                || !(confidenceLevel > 0 && confidenceLevel < 1)
                || !Double.isFinite(confidenceLevel))
            throw new IllegalArgumentException("invalid Cox optimization options");
    }

    public static CoxOptions defaults() {
        return new CoxOptions(CoxTies.EFRON, 50, 1e-9, 1e-7,
            20, 1e-10, 0.95);
    }

    public CoxOptions withTies(CoxTies value) {
        return new CoxOptions(value, maximumIterations, relativeTolerance,
            scoreTolerance, maximumStepHalvings, informationRidge,
            confidenceLevel);
    }
}
