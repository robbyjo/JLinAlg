/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

/** One feature's effect and method-specific empirical-Bayes inference. */
public record DifferentialResult(
        int featureIndex,
        double effect,
        double log2FoldChange,
        double standardError,
        double statistic,
        double degreesOfFreedom,
        double pValue,
        double baseMean,
        double rawVariance,
        double moderatedVariance,
        double dispersion,
        double meanPrecisionWeight,
        boolean converged,
        int iterations) { }
