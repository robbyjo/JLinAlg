/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.List;
import org.jlinalg.pipeline.VariantFilterResult;

/** Correlated-parametric-null SKAT-O omnibus result. */
public record SkatOResult(
        String setId,
        int requestedVariants,
        int includedVariants,
        List<Component> components,
        double minimumComponentPValue,
        double adjustedPValue,
        double log10AdjustedPValue,
        int simulations,
        long randomSeed,
        List<VariantFilterResult> excludedVariants) {
    public SkatOResult {
        components = List.copyOf(components);
        excludedVariants = List.copyOf(excludedVariants);
    }
    public double negativeLog10AdjustedPValue() {
        return -log10AdjustedPValue;
    }
    public record Component(double rho, SetTestResult result) { }
}
