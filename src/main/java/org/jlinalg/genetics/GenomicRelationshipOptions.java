/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

/** Variant filters used while constructing a genomic relationship matrix. */
public record GenomicRelationshipOptions(
        double minimumMinorAlleleFrequency,
        double minimumCallRate) {

    public GenomicRelationshipOptions {
        if (!Double.isFinite(minimumMinorAlleleFrequency)
                || minimumMinorAlleleFrequency < 0
                || minimumMinorAlleleFrequency >= 0.5)
            throw new IllegalArgumentException(
                "minimumMinorAlleleFrequency must lie in [0, 0.5)");
        if (!Double.isFinite(minimumCallRate)
                || minimumCallRate <= 0 || minimumCallRate > 1)
            throw new IllegalArgumentException(
                "minimumCallRate must lie in (0, 1]");
    }

    public static GenomicRelationshipOptions defaults() {
        return new GenomicRelationshipOptions(0.01, 0.95);
    }
}
