/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import org.jlinalg.pipeline.VariantRecord;

/** One explicitly oriented and weighted member of a variant set. */
public record WeightedVariant(
        VariantRecord variant, EffectAllele effectAllele, double weight) {
    public WeightedVariant {
        if (variant == null || effectAllele == null)
            throw new IllegalArgumentException("variant and effect allele are required");
        if (!(weight > 0) || !Double.isFinite(weight))
            throw new IllegalArgumentException(
                "variant-set weight must be finite and positive");
    }
}
