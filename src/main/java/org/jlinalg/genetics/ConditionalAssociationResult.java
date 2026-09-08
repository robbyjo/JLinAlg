/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import java.util.List;

/** SNP association inference; pValue is not a causal MR effect test. */
public record ConditionalAssociationResult(String variantId, double effect,
        double standardError, double statistic, double pValue,
        double confidenceLower, double confidenceUpper,
        double residualGenotypeVariance, List<String> conditionedOn) {
    public ConditionalAssociationResult {
        conditionedOn = List.copyOf(conditionedOn);
    }
}
