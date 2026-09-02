/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

final class DiscretePSplineBasisTest {
    @Test
    void uniqueRowsProduceExactExpandedCrossProducts() {
        double[] x = new double[120];
        double[] y = new double[120];
        double[] weights = new double[120];
        for (int row = 0; row < x.length; row++) {
            x[row] = (row % 8) / 7.0;
            y[row] = Math.sin(row);
            weights[row] = 0.5 + (row % 5) / 5.0;
        }
        PSplineTerm term = PSplineTerm.of("s(x)", x, 9);
        DiscretePSplineBasis discrete = DiscretePSplineBasis.compile(term);
        assertEquals(8, discrete.uniqueRows());
        assertArrayEquals(term.design(), discrete.expand(), 0.0);
        double[] expectedCross = new double[81];
        double[] expectedRight = new double[9];
        double[] basis = term.design();
        for (int row = 0; row < x.length; row++) {
            for (int first = 0; first < 9; first++) {
                expectedRight[first] += basis[row * 9 + first]
                    * weights[row] * y[row];
                for (int second = 0; second < 9; second++) {
                    expectedCross[first * 9 + second] +=
                        basis[row * 9 + first] * weights[row]
                            * basis[row * 9 + second];
                }
            }
        }
        assertArrayEquals(expectedCross, discrete.crossProduct(weights), 1e-12);
        assertArrayEquals(expectedRight,
            discrete.transposeMultiply(y, weights), 1e-12);
    }
}
