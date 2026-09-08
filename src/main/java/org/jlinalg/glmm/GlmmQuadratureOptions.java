/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Controls for scalar adaptive quadrature. The tolerance is absolute error
 * in the total log likelihood estimated by node refinements, not a rigorous bound. */
public record GlmmQuadratureOptions(int initialNodes, int maximumNodes,
        double quadratureTolerance, int maximumIterations, double gradientTolerance) {
    public GlmmQuadratureOptions {
        if (initialNodes < 1 || maximumNodes < initialNodes || maximumNodes > 512
                || !Double.isFinite(quadratureTolerance) || quadratureTolerance <= 0
                || maximumIterations < 1 || !Double.isFinite(gradientTolerance)
                || gradientTolerance <= 0)
            throw new IllegalArgumentException("invalid adaptive quadrature controls");
    }
    public static GlmmQuadratureOptions defaults() {
        return new GlmmQuadratureOptions(9, 257, 1e-10, 300, 1e-5);
    }
}
