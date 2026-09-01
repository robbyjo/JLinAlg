/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/**
 * Reusable Gaussian score projection for unrelated or related samples.
 * Variant rows are already aligned, oriented, imputed, and weighted.
 */
public interface GaussianSetTestNullModel {
    int observations();
    double degreesOfFreedom();
    SetTestScoreState score(double[][] variantRows);
}
