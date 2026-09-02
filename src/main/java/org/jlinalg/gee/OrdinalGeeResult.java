/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Arrays;
import java.util.Objects;

/** Cumulative-link ordinal GEE thresholds, slopes, and underlying fit. */
public final class OrdinalGeeResult {
    private final int categories;
    private final int predictors;
    private final boolean[] proportional;
    private final GeeResult fit;

    OrdinalGeeResult(int categories, GeeResult fit) {
        this(categories, fit.parameters() - categories + 1,
            allTrue(fit.parameters() - categories + 1), fit);
    }

    OrdinalGeeResult(
            int categories, int predictors,
            boolean[] proportional, GeeResult fit) {
        this.categories = categories;
        this.predictors = predictors;
        this.proportional = proportional.clone();
        this.fit = Objects.requireNonNull(fit, "fit");
    }

    public int categories() { return categories; }
    public double[] thresholds() {
        return Arrays.copyOf(fit.coefficients(), categories - 1);
    }
    public double[] coefficients() {
        double[] all = fit.coefficients();
        return Arrays.copyOfRange(all, categories - 1, all.length);
    }
    public int predictors() { return predictors; }
    public boolean[] proportional() { return proportional.clone(); }
    /** Shared slopes; non-proportional predictors are represented by NaN. */
    public double[] proportionalCoefficients() {
        double[] all = fit.coefficients();
        double[] result = new double[predictors];
        int source = categories - 1;
        for (int predictor = 0; predictor < predictors; predictor++) {
            if (proportional[predictor]) {
                result[predictor] = all[source++];
            } else {
                result[predictor] = Double.NaN;
                source += categories - 1;
            }
        }
        return result;
    }
    /** Row-major cutoff-by-predictor slope matrix. */
    public double[] cutoffSpecificCoefficients() {
        int cutoffs = categories - 1;
        double[] all = fit.coefficients();
        double[] result = new double[cutoffs * predictors];
        int source = cutoffs;
        for (int predictor = 0; predictor < predictors; predictor++) {
            if (proportional[predictor]) {
                for (int cutoff = 0; cutoff < cutoffs; cutoff++) {
                    result[cutoff * predictors + predictor] = all[source];
                }
                source++;
            } else {
                for (int cutoff = 0; cutoff < cutoffs; cutoff++) {
                    result[cutoff * predictors + predictor] = all[source++];
                }
            }
        }
        return result;
    }
    public GeeResult fit() { return fit; }

    private static boolean[] allTrue(int length) {
        boolean[] result = new boolean[length];
        Arrays.fill(result, true);
        return result;
    }
}
