/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.List;
import org.jlinalg.pipeline.VariantFilterResult;

/** Canonical six-component ACAT-O rare-variant omnibus result. */
public record AcatOResult(
        String setId,
        int requestedVariants,
        int includedVariants,
        double minorAlleleCountThreshold,
        List<Component> components,
        double statistic,
        double pValue,
        double log10PValue,
        List<VariantFilterResult> excludedVariants) {
    public AcatOResult {
        components = List.copyOf(components);
        excludedVariants = List.copyOf(excludedVariants);
    }
    public double negativeLog10PValue() { return -log10PValue; }
    public record Component(String name, double pValue) { }
}
