/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;

/** One exposure-outcome pair passing an xWAS MR significance screen. */
public record XwasMrHit(
        int exposureIndex,
        int outcomeIndex,
        String exposureId,
        String exposureLabel,
        String outcomeId,
        String outcomeLabel,
        String outcomeCategory,
        MrEstimate screeningEstimate,
        double negativeLog10PValue,
        MrAnalysisResult analysis,
        List<HarmonizedInstrument> harmonizedInstruments,
        List<HarmonizationExclusion> harmonizationExclusions) {
    /** Makes hit-only instruments and diagnostic exclusions immutable. */
    public XwasMrHit {
        harmonizedInstruments = List.copyOf(harmonizedInstruments);
        harmonizationExclusions = List.copyOf(harmonizationExclusions);
    }
}
