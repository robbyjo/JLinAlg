/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import java.util.List;

/** Warm-started sequence of penalized Gaussian regression fits. */
public final class PenalizedRegressionPath {
    private final List<PenalizedRegressionResult> fits;
    private final double[] lambdas;

    PenalizedRegressionPath(List<PenalizedRegressionResult> fits) {
        this.fits = List.copyOf(fits);
        this.lambdas = new double[fits.size()];
        for (int index = 0; index < fits.size(); index++) {
            lambdas[index] = fits.get(index).lambda();
        }
    }

    public List<PenalizedRegressionResult> fits() { return fits; }
    public double[] lambdas() { return lambdas.clone(); }
    public int size() { return fits.size(); }
    public PenalizedRegressionResult fit(int index) { return fits.get(index); }

    /** Returns the nearest fitted lambda on the log scale. */
    public PenalizedRegressionResult nearest(double lambda) {
        if (!Double.isFinite(lambda) || lambda < 0.0) {
            throw new IllegalArgumentException("lambda must be finite and nonnegative");
        }
        int best = 0;
        double bestDistance = distance(lambdas[0], lambda);
        for (int index = 1; index < lambdas.length; index++) {
            double candidate = distance(lambdas[index], lambda);
            if (candidate < bestDistance) {
                best = index;
                bestDistance = candidate;
            }
        }
        return fits.get(best);
    }

    private static double distance(double first, double second) {
        if (first == 0.0 || second == 0.0) {
            return Math.abs(first - second);
        }
        return Math.abs(Math.log(first) - Math.log(second));
    }
}
