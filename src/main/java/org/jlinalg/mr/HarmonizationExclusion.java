/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Objects;

/** One variant excluded during allele harmonization. */
public record HarmonizationExclusion(
        String variantId,
        HarmonizationExclusionReason reason,
        String detail) {
    /** Validates exclusion metadata. */
    public HarmonizationExclusion {
        Objects.requireNonNull(variantId, "variantId");
        Objects.requireNonNull(reason, "reason");
        Objects.requireNonNull(detail, "detail");
    }
}
