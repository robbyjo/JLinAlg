/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;

/** Instrument-level variance-explained directionality checks. */
public record SteigerResult(
        List<String> retainedVariants,
        List<String> reversedVariants,
        double exposureVarianceExplained,
        double outcomeVarianceExplained,
        boolean aggregateDirectionCorrect) {
    public SteigerResult {
        retainedVariants = List.copyOf(retainedVariants);
        reversedVariants = List.copyOf(reversedVariants);
    }
}
