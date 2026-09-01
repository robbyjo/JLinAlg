/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.List;
import org.jlinalg.pipeline.VariantFilterResult;

/** Common set-test result with membership audit and optional effect estimate. */
public record SetTestResult(
        String setId,
        String method,
        int requestedVariants,
        int includedVariants,
        double statistic,
        double beta,
        double standardError,
        double degreesOfFreedom,
        double pValue,
        double log10PValue,
        String pValueMethod,
        double[] eigenvalues,
        List<VariantFilterResult> excludedVariants) {
    public SetTestResult {
        eigenvalues = eigenvalues == null ? new double[0] : eigenvalues.clone();
        excludedVariants = List.copyOf(excludedVariants);
    }
    @Override public double[] eigenvalues() { return eigenvalues.clone(); }
    public double negativeLog10PValue() { return -log10PValue; }
}
