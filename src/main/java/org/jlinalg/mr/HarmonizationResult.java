/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;

/** Immutable retained instruments and structured harmonization exclusions. */
public record HarmonizationResult(
        List<HarmonizedInstrument> instruments,
        List<HarmonizationExclusion> exclusions) {
    /** Makes both result lists immutable. */
    public HarmonizationResult {
        instruments = List.copyOf(instruments);
        exclusions = List.copyOf(exclusions);
    }
}
