/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import jdistlib.T;

/**
 * Reusable Gaussian score projection for unrelated or related samples.
 * Variant rows are already aligned, oriented, imputed, and weighted.
 */
public interface GaussianSetTestNullModel extends SetTestScoreNullModel {
    double degreesOfFreedom();

    @Override
    default double burdenPValue(double statistic) {
        double logP = Math.log(2.0)
            + T.cumulative(-Math.abs(statistic), degreesOfFreedom(), true, true);
        return Math.min(1.0, Math.exp(logP));
    }

    @Override default String burdenPValueMethod() { return "student-t-score"; }
}
