/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;

/** The core univariable, uncorrelated-instrument summary-data MR analyses. */
public record MrAnalysisResult(
        List<WaldRatio> waldRatios,
        MrEstimate ivwFixed,
        MrEstimate ivwMultiplicativeRandom,
        MrEggerResult egger,
        MrEstimate weightedMedian,
        List<LeaveOneOutEstimate> leaveOneOut,
        double meanFStatistic,
        List<String> warnings) {
    /** Makes collection-valued results immutable. */
    public MrAnalysisResult {
        waldRatios = List.copyOf(waldRatios);
        leaveOneOut = List.copyOf(leaveOneOut);
        warnings = List.copyOf(warnings);
    }
}
