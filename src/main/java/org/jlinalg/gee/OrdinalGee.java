/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;

/** Proportional-odds ordinal GEE through cumulative binary indicators. */
public final class OrdinalGee {
    private OrdinalGee() { }

    /**
     * Fits {@code logit(P(Y <= k)) = threshold[k] - x beta}.
     * The covariate matrix must not include an intercept; category values are
     * zero based. Working independence is used so the sandwich covariance
     * accounts for both repeated observations and cumulative indicators.
     */
    public static OrdinalGeeResult fit(
            int[] response,
            double[][] covariates,
            int[] cluster,
            int[] repeated,
            int categories,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || covariates == null || cluster == null
                || options == null || backendPolicy == null
                || response.length != covariates.length
                || response.length != cluster.length
                || repeated != null && repeated.length != response.length
                || response.length == 0 || covariates[0] == null) {
            throw new IllegalArgumentException("ordinal GEE input dimensions are invalid");
        }
        if (categories < 3) {
            throw new IllegalArgumentException(
                "ordinal GEE requires at least three categories");
        }
        if (options.correlation() != GeeCorrelation.INDEPENDENCE) {
            throw new IllegalArgumentException(
                "ordinal GEE currently requires working independence");
        }
        int n = response.length;
        int predictors = covariates[0].length;
        int cutoffs = categories - 1;
        int augmentedRows = n * cutoffs;
        int columns = cutoffs + predictors;
        double[] binary = new double[augmentedRows];
        double[][] design = new double[augmentedRows][columns];
        int[] augmentedCluster = new int[augmentedRows];
        int[] augmentedRepeated = new int[augmentedRows];
        for (int row = 0; row < n; row++) {
            if (response[row] < 0 || response[row] >= categories) {
                throw new IllegalArgumentException(
                    "ordinal responses must be in [0, categories)");
            }
            if (covariates[row] == null || covariates[row].length != predictors) {
                throw new IllegalArgumentException("covariates must be rectangular");
            }
            int wave = repeated == null ? row : repeated[row];
            for (int cutoff = 0; cutoff < cutoffs; cutoff++) {
                int destination = row * cutoffs + cutoff;
                binary[destination] = response[row] <= cutoff ? 1.0 : 0.0;
                design[destination][cutoff] = 1.0;
                for (int column = 0; column < predictors; column++) {
                    design[destination][cutoffs + column] = -covariates[row][column];
                }
                augmentedCluster[destination] = cluster[row];
                augmentedRepeated[destination] = wave * cutoffs + cutoff;
            }
        }
        GeeResult fit = Gee.fit(binary, design, augmentedCluster,
            augmentedRepeated, GlmFamilies.binomial(), null, null,
            options, backendPolicy);
        return new OrdinalGeeResult(categories, fit);
    }
}
