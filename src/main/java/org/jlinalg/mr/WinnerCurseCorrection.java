/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import jdistlib.Normal;

/** Selection-adjusted normal maximum likelihood for a genome-wide hit. */
public final class WinnerCurseCorrection {
    private WinnerCurseCorrection() { }

    /** Corrects an effect selected because {@code |beta / se| >= threshold}. */
    public static double correct(double beta, double standardError, double threshold) {
        if (!(standardError > 0.0) || !(threshold > 0.0)
                || !Double.isFinite(beta) || !Double.isFinite(standardError)
                || !Double.isFinite(threshold)) {
            throw new IllegalArgumentException("effect, SE, and threshold are invalid");
        }
        double lower = Math.min(-Math.abs(beta) * 2.0, beta - 10.0 * standardError);
        double upper = Math.max(Math.abs(beta) * 2.0, beta + 10.0 * standardError);
        double best = beta;
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        for (int point = 0; point <= 20_000; point++) {
            double mean = lower + (upper - lower) * point / 20_000.0;
            double selected = Normal.cumulative(-threshold, mean / standardError,
                    1.0, true, false)
                + Normal.cumulative(threshold, mean / standardError,
                    1.0, false, false);
            double logLikelihood = Normal.density(beta, mean, standardError, true)
                - Math.log(Math.max(1e-300, selected));
            if (logLikelihood > bestLikelihood) {
                bestLikelihood = logLikelihood;
                best = mean;
            }
        }
        return best;
    }
}
