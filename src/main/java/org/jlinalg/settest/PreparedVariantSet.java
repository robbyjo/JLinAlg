/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.List;
import org.jlinalg.pipeline.VariantFilterResult;

/** Filtered, oriented, imputed variant-set data reusable across null models. */
public final class PreparedVariantSet {
    private final String id;
    private final int requestedVariants;
    private final double[][] dosages;
    private final double[] weights;
    private final List<VariantFilterResult> excludedVariants;

    PreparedVariantSet(
            String id, int requestedVariants, double[][] dosages,
            double[] weights, List<VariantFilterResult> excludedVariants) {
        this.id = id;
        this.requestedVariants = requestedVariants;
        this.dosages = dosages;
        this.weights = weights;
        this.excludedVariants = List.copyOf(excludedVariants);
    }

    public String id() { return id; }
    public int requestedVariants() { return requestedVariants; }
    public int includedVariants() { return dosages.length; }
    public int observations() { return dosages[0].length; }
    public List<VariantFilterResult> excludedVariants() {
        return excludedVariants;
    }

    double[][] dosagesView() { return dosages; }
    double[] weightsView() { return weights; }
}
