/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Objects;

/** One successfully screened exposure-outcome pair in deterministic grid order. */
public record XwasMrScreeningResult(
        int exposureIndex,
        int outcomeIndex,
        String exposureId,
        String exposureLabel,
        String outcomeId,
        String outcomeLabel,
        String outcomeCategory,
        MrEstimate estimate,
        double negativeLog10PValue,
        boolean thresholdPassed) {
    /** Validates pair metadata and the transformed screening p-value. */
    public XwasMrScreeningResult {
        if (exposureIndex < 0 || outcomeIndex < 0)
            throw new IllegalArgumentException("pair indices must be nonnegative");
        Objects.requireNonNull(exposureId, "exposureId");
        Objects.requireNonNull(exposureLabel, "exposureLabel");
        Objects.requireNonNull(outcomeId, "outcomeId");
        Objects.requireNonNull(outcomeLabel, "outcomeLabel");
        Objects.requireNonNull(outcomeCategory, "outcomeCategory");
        Objects.requireNonNull(estimate, "estimate");
        if (Double.isNaN(negativeLog10PValue)
                || negativeLog10PValue < 0.0)
            throw new IllegalArgumentException(
                "negativeLog10PValue must be nonnegative");
    }
}
