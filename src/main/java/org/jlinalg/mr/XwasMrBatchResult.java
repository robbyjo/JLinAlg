/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;

/** Deterministically ordered significant xWAS MR pairs and scan provenance. */
public record XwasMrBatchResult(
        List<XwasMrHit> hits,
        List<XwasMrFailure> failures,
        long totalPairs,
        long analyzablePairs,
        long belowThresholdPairs,
        long insufficientInstrumentPairs,
        int workersUsed,
        long elapsedNanoseconds) {
    /** Makes result collections immutable and validates accounting. */
    public XwasMrBatchResult {
        hits = List.copyOf(hits);
        failures = List.copyOf(failures);
        if (totalPairs < 0L || analyzablePairs < 0L
                || belowThresholdPairs < 0L
                || insufficientInstrumentPairs < 0L
                || workersUsed < 1 || elapsedNanoseconds < 0L)
            throw new IllegalArgumentException("invalid xWAS scan provenance");
    }
}
