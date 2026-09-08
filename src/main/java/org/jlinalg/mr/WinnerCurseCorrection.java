/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import jdistlib.Normal;

/** Selection-adjusted normal maximum likelihood for a genome-wide hit. */
public final class WinnerCurseCorrection {
    private WinnerCurseCorrection() { }

    /**
     * Corrects an effect selected because {@code |beta / se| >= threshold}.
     * Solves the strictly concave, two-sided truncated-normal likelihood score
     * in standardized units; an unselected observation is not a valid input.
     * Selection thresholds above 1e6 are rejected as outside numerical scope.
     */
    public static double correct(double beta, double standardError, double threshold) {
        if (!(standardError > 0.0) || !(threshold > 0.0)
                || !Double.isFinite(beta) || !Double.isFinite(standardError)
                || !Double.isFinite(threshold)) {
            throw new IllegalArgumentException("effect, SE, and threshold are invalid");
        }
        double observed = Math.abs(beta) / standardError;
        if (!Double.isFinite(observed) || observed < threshold || threshold > 1e6)
            throw new IllegalArgumentException("observation must satisfy selection and have a representable standardized effect");
        // Selection changes the estimate by less than double precision here.
        if (observed - threshold > 40) return beta;
        double lower = 0.0, upper = observed;
        // E[Z | |Z| >= threshold] increases strictly with the normal mean.
        // At mean zero it is zero; at mean observed it is >= observed.
        for (int iteration = 0; iteration < 160; iteration++) {
            double mean = lower + 0.5 * (upper - lower);
            if (mean == lower || mean == upper) break;
            double a = threshold - mean, b = threshold + mean;
            double logUpper = Normal.cumulative(a, 0, 1, false, true);
            double logLower = Normal.cumulative(b, 0, 1, false, true);
            double ratio = Math.exp(logLower - logUpper);
            double selectedMean = mean
                + (inverseMills(a) - ratio * inverseMills(b)) / (1 + ratio);
            if (selectedMean > observed) upper = mean;
            else lower = mean;
        }
        return Math.copySign((lower + 0.5 * (upper - lower)) * standardError, beta);
    }

    private static double inverseMills(double z) {
        if (z > 8) {
            // Laplace's continued fraction avoids subtracting enormous log
            // density and log survival values in the far normal tail.
            double fraction = 0;
            for (int k = 64; k >= 1; k--) fraction = k / (z + fraction);
            return z + fraction;
        }
        return Math.exp(Normal.density(z, 0, 1, true)
            - Normal.cumulative(z, 0, 1, false, true));
    }
}
