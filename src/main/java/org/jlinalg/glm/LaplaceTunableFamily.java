/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glm;

/**
 * Family likelihood parameters optimized with GLMM variance components.
 * Implementations are mutable fit state and must not be shared by concurrent
 * prepared fits.
 */
public interface LaplaceTunableFamily extends GlmFamily {
    /** Current unconstrained parameters, normally on a logarithmic scale. */
    double[] laplaceParameters();
    void setLaplaceParameters(double[] parameters);
    double minimumLaplaceParameter(int index);
    double maximumLaplaceParameter(int index);
}
