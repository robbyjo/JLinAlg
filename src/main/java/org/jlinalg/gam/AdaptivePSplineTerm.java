/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;

/** P-spline with baseline and spatially weighted difference penalties. */
public final class AdaptivePSplineTerm {
    private AdaptivePSplineTerm() { }

    /**
     * Creates two positive-semidefinite penalties. Their fitted combination
     * permits the roughness cost to vary monotonically across the coefficient index.
     */
    public static QuadraticSmoothTerm of(String name, PSplineTerm marginal) {
        int columns = marginal.basisDimension();
        int rows = columns - marginal.differenceOrder();
        double[] difference = marginal.penaltyFactor();
        double[] baseline = new double[columns * columns];
        double[] adaptive = new double[columns * columns];
        for (int first = 0; first < columns; first++) {
            for (int second = 0; second < columns; second++) {
                for (int row = 0; row < rows; row++) {
                    double product = difference[row * columns + first]
                        * difference[row * columns + second];
                    baseline[first * columns + second] += product;
                    adaptive[first * columns + second] +=
                        (row + 1.0) / (rows + 1.0) * product;
                }
            }
        }
        return new QuadraticSmoothTerm(name, marginal.observations(), columns,
            marginal.design(), List.of(baseline, adaptive));
    }
}
