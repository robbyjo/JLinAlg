/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Controls frequency-based resolution of palindromic SNPs. */
public record HarmonizationOptions(double alleleFrequencyTolerance) {
    /** Validates harmonization controls. */
    public HarmonizationOptions {
        if (!(alleleFrequencyTolerance > 0.0)
                || !(alleleFrequencyTolerance < 0.5)
                || !Double.isFinite(alleleFrequencyTolerance)) {
            throw new IllegalArgumentException(
                "alleleFrequencyTolerance must be finite and lie in (0, 0.5)");
        }
    }

    /** Returns defaults suitable for matching summary-data allele frequencies. */
    public static HarmonizationOptions defaults() {
        return new HarmonizationOptions(0.08);
    }
}
