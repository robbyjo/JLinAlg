/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import jdistlib.ChiSquare;
import org.jlinalg.reml.VarianceEstimation;

/** Likelihood-ratio comparisons for caller-verified nested ML models. */
public final class MixedModelComparison {
    private MixedModelComparison() { }

    public static MixedModelComparisonResult compare(
            LinearMixedModelResult reduced,
            LinearMixedModelResult full) {
        if (reduced == null || full == null)
            throw new IllegalArgumentException("both models are required");
        return compare(reduced.reml().logLikelihood(),
            reduced.beta().length + reduced.reml().varianceComponents().length,
            reduced.reml().varianceEstimation(),
            full.reml().logLikelihood(),
            full.beta().length + full.reml().varianceComponents().length,
            full.reml().varianceEstimation());
    }

    public static MixedModelComparisonResult compare(
            SparseLinearMixedModelResult reduced,
            SparseLinearMixedModelResult full) {
        if (reduced == null || full == null)
            throw new IllegalArgumentException("both models are required");
        return compare(reduced.logLikelihood(),
            reduced.beta().length + reduced.varianceComponents().length,
            reduced.varianceEstimation(), full.logLikelihood(),
            full.beta().length + full.varianceComponents().length,
            full.varianceEstimation());
    }

    private static MixedModelComparisonResult compare(
            double reducedLikelihood, int reducedParameters,
            VarianceEstimation reducedEstimation,
            double fullLikelihood, int fullParameters,
            VarianceEstimation fullEstimation) {
        if (reducedEstimation != VarianceEstimation.ML
                || fullEstimation != VarianceEstimation.ML)
            throw new IllegalArgumentException(
                "models with different fixed effects must be compared using ML");
        int degrees = fullParameters - reducedParameters;
        if (degrees < 1)
            throw new IllegalArgumentException(
                "full model must contain more estimated parameters");
        double statistic = Math.max(0.0,
            2.0 * (fullLikelihood - reducedLikelihood));
        double pValue = Math.min(1.0, ChiSquare.cumulative(
            statistic, degrees, false, false));
        return new MixedModelComparisonResult(reducedLikelihood,
            fullLikelihood, statistic, degrees, pValue);
    }
}
