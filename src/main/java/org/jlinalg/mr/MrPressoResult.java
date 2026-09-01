/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;

/** PRESSO-style global heterogeneity and outlier-corrected estimate. */
public record MrPressoResult(
        MrEstimate rawEstimate,
        MrEstimate outlierCorrectedEstimate,
        List<String> outlierVariants,
        double globalStatistic,
        double globalPValue,
        double distortionPercent) {
    public MrPressoResult { outlierVariants = List.copyOf(outlierVariants); }
}
