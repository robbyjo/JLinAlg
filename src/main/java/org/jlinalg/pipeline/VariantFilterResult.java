/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.List;

/** Variant statistics plus every applicable exclusion reason. */
public record VariantFilterResult(
        VariantRecord variant,
        VariantStatistics statistics,
        List<VariantFilterReason> reasons) {
    public VariantFilterResult {
        reasons = List.copyOf(reasons);
    }
    public boolean included() { return reasons.isEmpty(); }
}
