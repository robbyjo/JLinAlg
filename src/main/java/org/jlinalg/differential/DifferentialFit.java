/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

import java.util.List;

/** Complete differential-analysis result and fitted cross-feature prior. */
public record DifferentialFit(
        String method,
        List<DifferentialResult> results,
        double priorDegreesOfFreedom,
        double priorVariance,
        double[] sizeFactors) {
    public DifferentialFit {
        if (method == null || method.isBlank())
            throw new IllegalArgumentException("method is required");
        results = List.copyOf(results);
        sizeFactors = sizeFactors.clone();
    }
    @Override public double[] sizeFactors() { return sizeFactors.clone(); }
}
