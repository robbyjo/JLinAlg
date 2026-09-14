/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.List;

/** Estimated sample factors and feature loadings for an omics matrix. */
public record LatentFactorResult(
        double[][] factors,
        double[][] loadings,
        double[] varianceExplained,
        double[][] adjusted,
        int iterations,
        boolean converged,
        List<Double> objectiveHistory) {
    public LatentFactorResult {
        objectiveHistory = List.copyOf(objectiveHistory);
    }
}
