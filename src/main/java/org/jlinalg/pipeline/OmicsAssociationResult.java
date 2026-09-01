/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.List;

/** Ordered result of a bounded-memory TWAS, EWAS, or PWAS matrix scan. */
public record OmicsAssociationResult(
        long sourceFeatures,
        List<OmicsAssociationEstimate> estimates,
        List<AssociationPipelineFailure> failures) {
    public OmicsAssociationResult {
        if (sourceFeatures < 0)
            throw new IllegalArgumentException("source feature count cannot be negative");
        estimates = List.copyOf(estimates);
        failures = List.copyOf(failures);
    }
}
