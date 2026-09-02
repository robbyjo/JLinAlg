/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.RemlOptions;

/** Exact Gaussian GAM with inverse-variance weights and an optional offset. */
public final class WeightedGam {
    private WeightedGam() { }

    /** Fits by using diag(1 / weight) as the residual covariance basis. */
    public static GammResult fitGaussian(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            double[] weights,
            double[] offset,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || weights == null || weights.length != response.length) {
            throw new IllegalArgumentException("response and one weight per row are required");
        }
        double[] adjusted = MatrixOps.finiteCopy(response, "response");
        double[] offsets = offset == null ? new double[response.length]
            : MatrixOps.finiteCopy(offset, "offset");
        if (offsets.length != response.length) {
            throw new IllegalArgumentException("offset length must match response");
        }
        double[] residualCovariance = new double[response.length * response.length];
        for (int row = 0; row < response.length; row++) {
            if (!(weights[row] > 0.0) || !Double.isFinite(weights[row])) {
                throw new IllegalArgumentException("weights must be finite and positive");
            }
            adjusted[row] -= offsets[row];
            residualCovariance[row * response.length + row] = 1.0 / weights[row];
        }
        return Gamm.fitGaussian(adjusted, parametricDesign, smoothTerms,
            List.of(), residualCovariance, options, backendPolicy);
    }

    /** Fits with unit weights, useful for applying only an offset. */
    public static GammResult fitGaussianWithOffset(
            double[] response,
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            double[] offset) {
        double[] weights = new double[response.length];
        Arrays.fill(weights, 1.0);
        return fitGaussian(response, parametricDesign, smoothTerms,
            weights, offset, RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }
}
