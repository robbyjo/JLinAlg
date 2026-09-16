/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;

/** Parametric-bootstrap multivariate MR-PRESSO result. */
public record MultivariateMrPressoResult(
        MultivariateMrResult rawEstimate,
        MultivariateMrResult outlierCorrectedEstimate,
        double globalStatistic,
        double globalPValue,
        List<MultivariateMrOutlier> instrumentTests,
        List<String> outlierVariants,
        double[] distortionPercent,
        int simulations,
        long seed) {
    public MultivariateMrPressoResult {
        instrumentTests = List.copyOf(instrumentTests);
        outlierVariants = List.copyOf(outlierVariants);
        distortionPercent = distortionPercent.clone();
    }
    @Override public double[] distortionPercent() {
        return distortionPercent.clone();
    }
}
