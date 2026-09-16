/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.List;
import org.jlinalg.pipeline.VariantFilterResult;

/** Variant-level ACAT result, including any ultra-rare burden component. */
public record AcatVResult(
        String setId,
        int requestedVariants,
        int includedVariants,
        String weighting,
        double minorAlleleCountThreshold,
        int collapsedVariants,
        List<Component> components,
        double statistic,
        double pValue,
        double log10PValue,
        List<VariantFilterResult> excludedVariants) {
    public AcatVResult {
        components = List.copyOf(components);
        excludedVariants = List.copyOf(excludedVariants);
    }
    public double negativeLog10PValue() { return -log10PValue; }
    public record Component(
        String name, int variants, double pValue, double normalizedWeight) { }
}
