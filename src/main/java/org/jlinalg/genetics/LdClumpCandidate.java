/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

/** One variant considered for linkage-disequilibrium clumping. */
public record LdClumpCandidate(String variantId, double pValue, String group) {
    /** Validates the identifier, p-value, and group. */
    public LdClumpCandidate {
        if (variantId == null || variantId.isBlank())
            throw new IllegalArgumentException("variantId must not be blank");
        if (!Double.isFinite(pValue) || pValue < 0.0 || pValue > 1.0)
            throw new IllegalArgumentException("pValue must lie in [0, 1]");
        if (group == null || group.isBlank())
            throw new IllegalArgumentException("group must not be blank");
    }
}
