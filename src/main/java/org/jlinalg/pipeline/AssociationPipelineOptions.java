/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Bounded-memory controls shared by streaming association pipelines. */
public record AssociationPipelineOptions(
        int variantBlockSize,
        VariantFilterOptions variantFilter) {
    public AssociationPipelineOptions {
        if (variantBlockSize < 1)
            throw new IllegalArgumentException("variant block size must be positive");
        if (variantFilter == null)
            throw new IllegalArgumentException("variant filter is required");
    }

    public static AssociationPipelineOptions defaults() {
        return new AssociationPipelineOptions(
            1024, VariantFilterOptions.defaults());
    }
}
