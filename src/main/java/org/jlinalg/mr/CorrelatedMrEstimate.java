/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;

/** A generalized MR estimate and the linear-algebra backend that produced it. */
public record CorrelatedMrEstimate(
        MrEstimate estimate,
        BackendProvenance backend) {
    /** Validates result members. */
    public CorrelatedMrEstimate {
        Objects.requireNonNull(estimate, "estimate");
        Objects.requireNonNull(backend, "backend");
    }
}
