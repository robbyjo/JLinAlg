/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

/**
 * Matrix-free finite-denominator-DF calculations for sparse REML callers.
 * The derivative vector is supplied by the sparse factorization layer, so the
 * final contraction never materializes observation-scale covariance matrices.
 */
public final class SparseFiniteDf {
    private SparseFiniteDf() { }

    /** Satterthwaite DF from a scalar contrast variance and its derivatives. */
    public static double satterthwaite(double variance,
            double[] varianceDerivatives, double[] varianceParameterCovariance) {
        if (!(variance > 0.0) || !Double.isFinite(variance)
                || varianceDerivatives == null || varianceParameterCovariance == null
                || varianceParameterCovariance.length != varianceDerivatives.length
                    * varianceDerivatives.length)
            throw new IllegalArgumentException("finite-DF inputs are invalid");
        double uncertainty = 0.0;
        int n = varianceDerivatives.length;
        for (int row = 0; row < n; row++) for (int column = 0; column < n; column++)
            uncertainty += varianceDerivatives[row] * varianceParameterCovariance[row * n + column]
                * varianceDerivatives[column];
        if (!(uncertainty > 0.0) || !Double.isFinite(uncertainty))
            return Double.POSITIVE_INFINITY;
        return Math.max(1.0, 2.0 * variance * variance / uncertainty);
    }

    /** Kenward-Roger-style covariance inflation from an additive adjustment. */
    public static double adjustedVariance(double modelVariance,
            double adjustment) {
        if (!(modelVariance > 0.0) || !Double.isFinite(modelVariance)
                || adjustment < 0.0 || !Double.isFinite(adjustment))
            throw new IllegalArgumentException("variance adjustment is invalid");
        return modelVariance + adjustment;
    }
}
