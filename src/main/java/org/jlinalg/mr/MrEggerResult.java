/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** MR-Egger causal slope, pleiotropy intercept, and instrument-strength diagnostic. */
public record MrEggerResult(
        MrEstimate slope,
        double intercept,
        double interceptStandardError,
        double interceptStatistic,
        double interceptPValue,
        double interceptConfidenceLower,
        double interceptConfidenceUpper,
        double iSquaredGx) {
}
