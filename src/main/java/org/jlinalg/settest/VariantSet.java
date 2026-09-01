/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.HashSet;
import java.util.List;

/** Named gene, region, pathway, or caller-defined variant set. */
public record VariantSet(String id, List<WeightedVariant> variants) {
    public VariantSet {
        if (id == null || id.isBlank())
            throw new IllegalArgumentException("variant-set ID is required");
        id = id.trim();
        variants = List.copyOf(variants);
        if (variants.isEmpty())
            throw new IllegalArgumentException("variant set cannot be empty");
        HashSet<String> seen = new HashSet<>();
        for (WeightedVariant variant : variants)
            if (!seen.add(variant.variant().id()))
                throw new IllegalArgumentException(
                    "variant set contains a duplicate ID: "
                    + variant.variant().id());
    }
}
