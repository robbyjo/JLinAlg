/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

/** PLINK/TwoSampleMR-compatible LD-clumping thresholds. */
public record LdClumpOptions(
        int windowKilobases,
        double rSquaredThreshold,
        double indexPValueThreshold) {
    /** TwoSampleMR {@code clump_data()} defaults. */
    public static LdClumpOptions defaults() {
        return new LdClumpOptions(10_000, 0.001, 1.0);
    }

    /** Validates all thresholds. */
    public LdClumpOptions {
        if (windowKilobases < 1)
            throw new IllegalArgumentException(
                "windowKilobases must be positive");
        if (!Double.isFinite(rSquaredThreshold)
                || rSquaredThreshold < 0.0 || rSquaredThreshold > 1.0)
            throw new IllegalArgumentException(
                "rSquaredThreshold must lie in [0, 1]");
        if (!Double.isFinite(indexPValueThreshold)
                || indexPValueThreshold < 0.0 || indexPValueThreshold > 1.0)
            throw new IllegalArgumentException(
                "indexPValueThreshold must lie in [0, 1]");
    }
}
