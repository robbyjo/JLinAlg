/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.Test;

final class AdditionalSmoothBasesTest {
    @Test
    void cyclicBasisAgreesAtPeriodBoundary() {
        QuadraticSmoothTerm term = CyclicSplineTerm.of(
            "s(hour,bs=cc)", new double[] {0.0, 6.0, 12.0, 18.0, 24.0},
            3, 24.0);
        double[] design = term.design();
        for (int column = 0; column < term.columns(); column++) {
            assertEquals(design[column], design[4 * term.columns() + column], 1e-12);
        }
    }

    @Test
    void thinPlateRandomEffectAndMrfExposeQuadraticPenalties() {
        int observations = 25;
        double[] x = new double[observations];
        double[] z = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = row % 5;
            z[row] = row / 5;
        }
        QuadraticSmoothTerm thin = ThinPlateSplineTerm.of2d(
            "s(x,z,bs=tp)", x, z, 7);
        assertEquals(10, thin.columns());
        QuadraticSmoothTerm random = RandomEffectSmoothTerm.of(
            "s(group,bs=re)", List.of("a", "b", "a", "c"));
        assertEquals(3, random.columns());
        int[] nodes = {0, 1, 2, 0, 1, 2};
        QuadraticSmoothTerm mrf = MarkovRandomFieldSmoothTerm.of("region",
            nodes, 3, List.of(new MarkovRandomFieldSmoothTerm.Edge(0, 1, 1.0),
                new MarkovRandomFieldSmoothTerm.Edge(1, 2, 1.0)));
        double[] penalty = mrf.penalties().get(0);
        for (int row = 0; row < 3; row++) {
            double sum = 0.0;
            for (int column = 0; column < 3; column++) {
                sum += penalty[row * 3 + column];
            }
            assertEquals(0.0, sum, 1e-12);
        }
        assertTrue(thin.penaltyCount() == 1 && mrf.penaltyCount() == 1);
    }
}
