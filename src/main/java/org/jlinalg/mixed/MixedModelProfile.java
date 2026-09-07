/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.function.DoubleUnaryOperator;
import jdistlib.ChiSquare;

/** Generic deterministic profile-likelihood interval solver for mixed-model refits. */
public final class MixedModelProfile {
    private MixedModelProfile() { }

    /**
     * Profiles a scalar variance or fixed-effect parameter. The callback must
     * return the maximized log likelihood with the parameter held fixed.
     */
    public static ProfileLikelihoodInterval interval(
            DoubleUnaryOperator profiledLogLikelihood, double estimate,
            double maximumLogLikelihood, double confidenceLevel,
            double lowerBound, double upperBound, int gridSize) {
        if (profiledLogLikelihood == null || !Double.isFinite(estimate)
                || !Double.isFinite(maximumLogLikelihood) || !(confidenceLevel > 0 && confidenceLevel < 1)
                || !(lowerBound < estimate) || !(upperBound > estimate) || gridSize < 8)
            throw new IllegalArgumentException("invalid profile-likelihood interval inputs");
        double cutoff = 0.5 * ChiSquare.quantile(confidenceLevel, 1.0, true, false);
        double target = maximumLogLikelihood - cutoff;
        double lower = Double.NaN, upper = Double.NaN;
        double previousX = estimate, previousValue = maximumLogLikelihood;
        for (int i = 1; i <= gridSize; i++) {
            double x = estimate - (estimate - lowerBound) * i / gridSize;
            double value = profiledLogLikelihood.applyAsDouble(x);
            if (Double.isFinite(value) && value < target) { lower = bisect(profiledLogLikelihood, previousX, x, target); break; }
            previousX = x; previousValue = value;
        }
        previousX = estimate; previousValue = maximumLogLikelihood;
        for (int i = 1; i <= gridSize; i++) {
            double x = estimate + (upperBound - estimate) * i / gridSize;
            double value = profiledLogLikelihood.applyAsDouble(x);
            if (Double.isFinite(value) && value < target) { upper = bisect(profiledLogLikelihood, previousX, x, target); break; }
            previousX = x; previousValue = value;
        }
        return new ProfileLikelihoodInterval(estimate, lower, upper, cutoff,
            Double.isFinite(lower), Double.isFinite(upper));
    }

    private static double bisect(DoubleUnaryOperator f, double inside,
            double outside, double target) {
        double left = inside, right = outside;
        for (int i = 0; i < 80; i++) {
            double middle = 0.5 * (left + right);
            double value = f.applyAsDouble(middle);
            if (!Double.isFinite(value) || value < target) right = middle; else left = middle;
        }
        return 0.5 * (left + right);
    }
}
