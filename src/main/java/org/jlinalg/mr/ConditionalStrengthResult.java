/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;

/** Cross-exposure-covariance-aware conditional instrument-strength diagnostics. */
public record ConditionalStrengthResult(List<String> exposureNames,
        double[] fStatistics, double[] qStatistics, int degreesOfFreedom) {
    public ConditionalStrengthResult {
        exposureNames = List.copyOf(exposureNames);
        fStatistics = fStatistics.clone();
        qStatistics = qStatistics.clone();
        if (exposureNames.size() != fStatistics.length
                || qStatistics.length != fStatistics.length
                || degreesOfFreedom < 1)
            throw new IllegalArgumentException("conditional-strength result dimensions are invalid");
    }
    @Override public double[] fStatistics() { return fStatistics.clone(); }
    @Override public double[] qStatistics() { return qStatistics.clone(); }
}
