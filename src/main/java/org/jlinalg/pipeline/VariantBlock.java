/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.List;

/** Ordered block of variants read from one source. */
public record VariantBlock(long firstSourceIndex, List<VariantRecord> variants) {
    public VariantBlock {
        if (firstSourceIndex < 0)
            throw new IllegalArgumentException("first source index cannot be negative");
        variants = List.copyOf(variants);
        if (variants.isEmpty())
            throw new IllegalArgumentException("variant block cannot be empty");
    }
}
