/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import org.jlinalg.internal.MatrixOps;

/** Low-rank two-dimensional thin-plate radial basis smooth. */
public final class ThinPlateSplineTerm {
    private ThinPlateSplineTerm() { }

    /** Creates affine null-space columns plus evenly sampled radial knots. */
    public static QuadraticSmoothTerm of2d(
            String name, double[] first, double[] second, int radialKnots) {
        double[] x = MatrixOps.finiteCopy(first, "first covariate");
        double[] z = MatrixOps.finiteCopy(second, "second covariate");
        if (x.length != z.length || x.length < 4 || radialKnots < 1
                || radialKnots > x.length) {
            throw new IllegalArgumentException("invalid thin-plate dimensions");
        }
        int columns = 3 + radialKnots;
        double[] design = new double[x.length * columns];
        double[] penalty = new double[columns * columns];
        for (int row = 0; row < x.length; row++) {
            design[row * columns] = 1.0;
            design[row * columns + 1] = x[row];
            design[row * columns + 2] = z[row];
            for (int knot = 0; knot < radialKnots; knot++) {
                int source = radialKnots == 1 ? x.length / 2
                    : (int) Math.round(knot * (x.length - 1.0) / (radialKnots - 1.0));
                double dx = x[row] - x[source];
                double dz = z[row] - z[source];
                double radiusSquared = dx * dx + dz * dz;
                design[row * columns + 3 + knot] = radiusSquared == 0.0
                    ? 0.0 : 0.5 * radiusSquared * Math.log(radiusSquared);
            }
        }
        for (int column = 3; column < columns; column++) {
            penalty[column * columns + column] = 1.0;
        }
        return new QuadraticSmoothTerm(name, x.length, columns,
            design, List.of(penalty));
    }
}
