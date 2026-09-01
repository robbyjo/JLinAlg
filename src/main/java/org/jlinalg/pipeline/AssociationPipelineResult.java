/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.List;

/** Ordered tested results and complete QC/failure audit from a streamed scan. */
public record AssociationPipelineResult(
        long sourceVariants,
        List<AssociationPipelineEstimate> estimates,
        List<VariantFilterResult> excludedVariants,
        List<AssociationPipelineFailure> failures) {
    public AssociationPipelineResult {
        if (sourceVariants < 0)
            throw new IllegalArgumentException("source variant count cannot be negative");
        estimates = List.copyOf(estimates);
        excludedVariants = List.copyOf(excludedVariants);
        failures = List.copyOf(failures);
    }
    public long testedVariants() { return estimates.size(); }
    public long excludedVariantCount() { return excludedVariants.size(); }
}
