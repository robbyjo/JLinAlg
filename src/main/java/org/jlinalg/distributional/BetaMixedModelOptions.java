/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import org.jlinalg.glmm.GlmmLaplaceOptions;

/** Controls for sparse beta mixed-model fitting. */
public record BetaMixedModelOptions(
        GlmmLaplaceOptions laplace,
        double initialPrecision,
        double minimumPrecision,
        double maximumPrecision) {
    public BetaMixedModelOptions {
        if (laplace == null
                || !(initialPrecision > 0.0) || !Double.isFinite(initialPrecision)
                || !(minimumPrecision > 0.0)
                || !(maximumPrecision > minimumPrecision)
                || initialPrecision < minimumPrecision
                || initialPrecision > maximumPrecision
                || !Double.isFinite(maximumPrecision)) {
            throw new IllegalArgumentException("invalid beta mixed-model controls");
        }
    }

    public static BetaMixedModelOptions defaults() {
        return new BetaMixedModelOptions(
            GlmmLaplaceOptions.defaults(), 10.0, 1e-6, 1e8);
    }
}
