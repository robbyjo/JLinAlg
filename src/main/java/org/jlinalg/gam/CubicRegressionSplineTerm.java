/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.internal.MatrixOps;

/** Low-rank cubic radial regression spline with an affine null space. */
public final class CubicRegressionSplineTerm {
    private CubicRegressionSplineTerm() { }

    /** Creates affine columns plus cubic radial basis functions at quantile knots. */
    public static QuadraticSmoothTerm of(
            String name, double[] covariate, int radialKnots) {
        double[] values = MatrixOps.finiteCopy(covariate, "covariate");
        if (values.length < 4 || radialKnots < 1 || radialKnots > values.length) {
            throw new IllegalArgumentException("invalid cubic regression spline dimensions");
        }
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        int columns = 2 + radialKnots;
        double[] design = new double[values.length * columns];
        double[] penalty = new double[columns * columns];
        for (int row = 0; row < values.length; row++) {
            design[row * columns] = 1.0;
            design[row * columns + 1] = values[row];
            for (int knot = 0; knot < radialKnots; knot++) {
                int index = radialKnots == 1 ? sorted.length / 2
                    : (int) Math.round(knot * (sorted.length - 1.0)
                        / (radialKnots - 1.0));
                design[row * columns + 2 + knot] =
                    Math.pow(Math.abs(values[row] - sorted[index]), 3.0);
            }
        }
        for (int column = 2; column < columns; column++) {
            penalty[column * columns + column] = 1.0;
        }
        return new QuadraticSmoothTerm(name, values.length, columns,
            design, List.of(penalty));
    }
}
