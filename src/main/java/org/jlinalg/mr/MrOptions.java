/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Immutable controls for a summary-data MR analysis. */
public record MrOptions(
        double confidenceLevel,
        int weightedMedianBootstrapReplicates,
        long randomSeed) {
    /** Validates analysis controls. */
    public MrOptions {
        if (!(confidenceLevel > 0.0) || !(confidenceLevel < 1.0)
                || !Double.isFinite(confidenceLevel)) {
            throw new IllegalArgumentException(
                "confidenceLevel must be finite and lie in (0, 1)");
        }
        if (weightedMedianBootstrapReplicates < 2) {
            throw new IllegalArgumentException(
                "weightedMedianBootstrapReplicates must be at least 2");
        }
    }

    /** Returns reproducible defaults. */
    public static MrOptions defaults() {
        return new MrOptions(0.95, 1000, 20260831L);
    }
}
