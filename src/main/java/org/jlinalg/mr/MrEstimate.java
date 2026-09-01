/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** One causal estimate with first-order inference and heterogeneity diagnostics. */
public record MrEstimate(
        MrMethod method,
        double estimate,
        double standardError,
        double statistic,
        double pValue,
        double confidenceLower,
        double confidenceUpper,
        double cochranQ,
        int heterogeneityDegreesOfFreedom,
        double heterogeneityPValue,
        double iSquared,
        double residualDispersion,
        int instrumentCount) {
}
