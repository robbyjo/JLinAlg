/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

/** Marginal association on the allele orientation of the supplied signed LD. */
public record ConditionalAssociation(String variantId, double marginalEffect,
        double marginalStandardError, String group) {
    public ConditionalAssociation {
        if (variantId == null || variantId.isBlank() || group == null || group.isBlank()
                || !Double.isFinite(marginalEffect) || !Double.isFinite(marginalStandardError)
                || !(marginalStandardError > 0.0)
                || !Double.isFinite(marginalEffect / marginalStandardError))
            throw new IllegalArgumentException("finite marginal effects, positive SEs and identifiers are required");
    }
}
