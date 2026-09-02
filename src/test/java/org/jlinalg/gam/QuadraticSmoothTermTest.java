/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

final class QuadraticSmoothTermTest {
    @Test
    void tensorProductHasTwoSymmetricMarginalPenalties() {
        int observations = 72;
        double[] first = new double[observations];
        double[] second = new double[observations];
        for (int row = 0; row < observations; row++) {
            first[row] = (row % 9) / 8.0;
            second[row] = (row / 9) / 7.0;
        }
        QuadraticSmoothTerm term = TensorProductPSplineTerm.of(
            "te(x,z)", first, second, 7, 6);
        assertEquals(observations, term.observations());
        assertEquals(42, term.columns());
        assertEquals(2, term.penaltyCount());
        for (double[] penalty : term.penalties()) {
            for (int row = 0; row < term.columns(); row++) {
                assertTrue(penalty[row * term.columns() + row] >= 0.0);
                for (int column = 0; column < row; column++) {
                    assertEquals(penalty[row * term.columns() + column],
                        penalty[column * term.columns() + row], 1e-12);
                }
            }
        }
    }

    @Test
    void arbitraryQuadraticPenaltyCompilesToDiagonalForm() {
        int observations = 90;
        double[] first = new double[observations];
        double[] second = new double[observations];
        for (int row = 0; row < observations; row++) {
            first[row] = ((row * 17) % observations) / 89.0;
            second[row] = ((row * 31) % observations) / 89.0;
        }
        QuadraticSmoothTerm term = TensorProductPSplineTerm.of(
            "surface", first, second, 6, 6);
        PenalizedPredictor predictor = QuadraticPenalizedPredictor.compile(
            intercept(observations), List.of(term),
            List.of(new double[] {2.0, 7.0}), BackendPolicy.CPU);
        assertEquals(observations, predictor.observations());
        assertEquals(term.columns(), predictor.columns());
        int positive = 0;
        for (double value : predictor.penaltyDiagonal()) {
            if (value > 0.0) positive++;
        }
        assertTrue(positive > 0);
        assertTrue(positive < predictor.columns());
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
