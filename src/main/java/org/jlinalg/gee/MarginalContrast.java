/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Link-scale marginal contrast and optional exponentiated effect. */
public record MarginalContrast(
        double estimate,
        double standardError,
        double statistic,
        double pValue,
        double confidenceLower,
        double confidenceUpper,
        double exponentiatedEstimate,
        double exponentiatedLower,
        double exponentiatedUpper) {
}
