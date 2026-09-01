/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;

/** A generalized MR-Egger result and its linear-algebra backend provenance. */
public record CorrelatedMrEggerResult(
        MrEggerResult estimate,
        BackendProvenance backend) {
    /** Validates result members. */
    public CorrelatedMrEggerResult {
        Objects.requireNonNull(estimate, "estimate");
        Objects.requireNonNull(backend, "backend");
    }
}
