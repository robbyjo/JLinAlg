/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

/** Reproducible multiple-imputation chain controls. */
public record MiceOptions(int imputations, int iterations,
        int predictiveMeanDonors, long seed, double ridge) {
    public MiceOptions {
        if (imputations < 2) throw new IllegalArgumentException("at least two imputations are required");
        if (iterations < 1) throw new IllegalArgumentException("iterations must be positive");
        if (predictiveMeanDonors < 1)
            throw new IllegalArgumentException("predictiveMeanDonors must be positive");
        if (!(ridge > 0.0) || !Double.isFinite(ridge))
            throw new IllegalArgumentException("ridge must be finite and positive");
    }
    public static MiceOptions defaults() { return new MiceOptions(5, 10, 5, 1L, 1e-6); }
}
