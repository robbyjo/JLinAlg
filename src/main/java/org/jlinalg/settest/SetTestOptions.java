/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import org.jlinalg.pipeline.VariantFilterOptions;

/** Variant membership, missingness, and SKAT-O calibration controls. */
public record SetTestOptions(
        VariantFilterOptions variantFilter,
        SetTestMissingPolicy missingPolicy,
        double[] skatORhoGrid,
        int skatOSimulations,
        long randomSeed) {
    public SetTestOptions {
        if (variantFilter == null || missingPolicy == null)
            throw new IllegalArgumentException(
                "variant filter and missing policy are required");
        if (skatORhoGrid == null || skatORhoGrid.length < 2)
            throw new IllegalArgumentException(
                "SKAT-O rho grid requires at least two values");
        skatORhoGrid = skatORhoGrid.clone();
        if (skatORhoGrid[0] != 0
                || skatORhoGrid[skatORhoGrid.length - 1] != 1)
            throw new IllegalArgumentException(
                "SKAT-O rho grid must begin at zero and end at one");
        for (int index = 0; index < skatORhoGrid.length; index++) {
            double rho = skatORhoGrid[index];
            if (!Double.isFinite(rho) || rho < 0 || rho > 1
                    || (index > 0 && rho <= skatORhoGrid[index - 1]))
                throw new IllegalArgumentException(
                    "SKAT-O rho values must increase strictly within [0,1]");
        }
        if (skatOSimulations < 1)
            throw new IllegalArgumentException(
                "SKAT-O simulations must be positive");
    }
    @Override public double[] skatORhoGrid() { return skatORhoGrid.clone(); }

    public static SetTestOptions defaults() {
        return new SetTestOptions(VariantFilterOptions.defaults(),
            SetTestMissingPolicy.MEAN_IMPUTE,
            new double[] {0, 0.25, 0.5, 0.75, 1}, 10_000, 20260901L);
    }
}
