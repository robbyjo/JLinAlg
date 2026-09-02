/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import org.jlinalg.internal.MatrixOps;

/** Public stationary AR/MA/ARMA residual-correlation builders for mixed models. */
public final class ArmaCorrelation {
    private ArmaCorrelation() { }

    /** Builds a stationary correlation matrix from conventional AR and MA coefficients. */
    public static double[] of(
            int observations,
            double[] autoregressive,
            double[] movingAverage) {
        double[] ar = autoregressive == null ? new double[0]
            : MatrixOps.finiteCopy(autoregressive, "autoregressive");
        double[] ma = movingAverage == null ? new double[0]
            : MatrixOps.finiteCopy(movingAverage, "movingAverage");
        validateStationary(ar);
        return ArimaMath.correlationMatrix(observations, ar, ma);
    }

    /** Convenience AR(1) correlation with entry phi^|i-j|. */
    public static double[] ar1(int observations, double phi) {
        if (!Double.isFinite(phi) || Math.abs(phi) >= 1.0) {
            throw new IllegalArgumentException("AR(1) phi must lie strictly in (-1,1)");
        }
        return of(observations, new double[] {phi}, null);
    }

    private static void validateStationary(double[] ar) {
        if (ar.length == 1 && Math.abs(ar[0]) >= 1.0) {
            throw new IllegalArgumentException("AR(1) coefficient is not stationary");
        }
        // The covariance expansion detects divergent higher-order polynomials
        // without rejecting stable AR(p) models whose absolute sum exceeds one.
    }
}
