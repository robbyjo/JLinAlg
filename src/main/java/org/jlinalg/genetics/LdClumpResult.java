/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import java.util.List;

/** Retained index variants and structured clumping exclusions. */
public record LdClumpResult(
        List<LdClumpCandidate> retained,
        List<LdClumpExclusion> exclusions) {
    /** Defensively copies result lists. */
    public LdClumpResult {
        retained = List.copyOf(retained);
        exclusions = List.copyOf(exclusions);
    }
}
