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
                || !Double.isFinite(lowerBound) || !Double.isFinite(upperBound)
                || !(lowerBound <= estimate) || !(upperBound >= estimate)
                || !(lowerBound < upperBound) || gridSize < 8)
            throw new IllegalArgumentException("invalid profile-likelihood interval inputs");
        return interval(profiledLogLikelihood,estimate,maximumLogLikelihood,
            confidenceLevel,lowerBound,upperBound,gridSize,
            0.5 * ChiSquare.quantile(confidenceLevel, 1.0, true, false));
    }

    /**
     * One-sided variance-boundary profile using the 50:50 point-mass/chi-square
     * likelihood-ratio law. The confidence level must exceed one half.
     */
    public static ProfileLikelihoodInterval boundaryInterval(
            DoubleUnaryOperator profiledLogLikelihood, double estimate,
            double maximumLogLikelihood, double confidenceLevel,
            double lowerBound, double upperBound, int gridSize) {
        if (!(confidenceLevel > .5))
            throw new IllegalArgumentException("boundary profile confidence must exceed one half");
        return interval(profiledLogLikelihood,estimate,maximumLogLikelihood,
            confidenceLevel,lowerBound,upperBound,gridSize,
            0.5 * ChiSquare.quantile(2*confidenceLevel-1,1,true,false));
    }

    private static ProfileLikelihoodInterval interval(
            DoubleUnaryOperator profiledLogLikelihood, double estimate,
            double maximumLogLikelihood, double confidenceLevel,
            double lowerBound, double upperBound, int gridSize,
            double cutoff) {
        double target = maximumLogLikelihood - cutoff;
        double lower = estimate == lowerBound ? lowerBound : Double.NaN;
        double upper = estimate == upperBound ? upperBound : Double.NaN;
        boolean lowerCrossing = false, upperCrossing = false;
        double previousX = estimate;
        for (int i = 1; lowerBound < estimate && i <= gridSize; i++) {
            double x = estimate - (estimate - lowerBound) * i / gridSize;
            double value = profiledLogLikelihood.applyAsDouble(x);
            requireFinite(value);
            if (value < target) { lower = bisect(profiledLogLikelihood, previousX, x, target); lowerCrossing = true; break; }
            previousX = x;
        }
        previousX = estimate;
        for (int i = 1; upperBound > estimate && i <= gridSize; i++) {
            double x = estimate + (upperBound - estimate) * i / gridSize;
            double value = profiledLogLikelihood.applyAsDouble(x);
            requireFinite(value);
            if (value < target) { upper = bisect(profiledLogLikelihood, previousX, x, target); upperCrossing = true; break; }
            previousX = x;
        }
        return new ProfileLikelihoodInterval(estimate, lower, upper, cutoff,
            lowerCrossing, upperCrossing);
    }

    private static double bisect(DoubleUnaryOperator f, double inside,
            double outside, double target) {
        double left = inside, right = outside;
        for (int i = 0; i < 80; i++) {
            double middle = 0.5 * (left + right);
            if (Math.abs(left - right) < 1e-8 * Math.max(1, Math.abs(middle))) break;
            double value = f.applyAsDouble(middle);
            requireFinite(value);
            if (value < target) right = middle; else left = middle;
        }
        return 0.5 * (left + right);
    }

    private static void requireFinite(double value) {
        if (!Double.isFinite(value))
            throw new IllegalStateException("profile refit returned a nonfinite likelihood");
    }
}
