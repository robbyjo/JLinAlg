/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Constant-memory counts returned by an incremental association scan. */
public record AssociationPipelineSummary(
        long sourceVariants,
        long testedVariants,
        long excludedVariants,
        long failures) {
    public AssociationPipelineSummary {
        if (sourceVariants < 0 || testedVariants < 0
                || excludedVariants < 0 || failures < 0)
            throw new IllegalArgumentException("pipeline counts cannot be negative");
    }
}
