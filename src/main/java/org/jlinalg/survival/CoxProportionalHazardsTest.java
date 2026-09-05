/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Log-time score tests for coefficient-specific and global PH departures. */
public record CoxProportionalHazardsTest(
        double[] statistics,
        double[] pValues,
        double globalStatistic,
        int globalDegreesOfFreedom,
        double globalPValue,
        int events,
        String timeTransform) {
    public CoxProportionalHazardsTest {
        statistics = statistics.clone();
        pValues = pValues.clone();
    }

    @Override public double[] statistics() { return statistics.clone(); }
    @Override public double[] pValues() { return pValues.clone(); }
}
