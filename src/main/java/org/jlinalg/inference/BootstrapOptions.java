/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

/** Deterministic parametric-bootstrap controls. */
public record BootstrapOptions(
        int simulations,
        double confidenceLevel,
        long randomSeed,
        int parallelism) {
    public BootstrapOptions {
        if (simulations < 2)
            throw new IllegalArgumentException(
                "at least two bootstrap simulations are required");
        if (!(confidenceLevel > 0 && confidenceLevel < 1)
                || !Double.isFinite(confidenceLevel))
            throw new IllegalArgumentException(
                "confidence level must be finite and in (0,1)");
        if (parallelism < 1)
            throw new IllegalArgumentException("parallelism must be positive");
    }

    /** Accelerator-friendly defaults; raise parallelism deliberately for CPU. */
    public static BootstrapOptions defaults() {
        return new BootstrapOptions(999, 0.95, 20260901L, 1);
    }
}
