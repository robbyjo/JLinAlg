/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

final class MultipleImputationTest {
    @Test void mixedChainsAreReproducibleAndPreserveObservedCells() {
        double missing = Double.NaN;
        double[][] data = {
            {20, 0, 1}, {21, 0, 1}, {22, 0, 2}, {missing, 1, 2},
            {30, 1, missing}, {31, 1, 2}, {32, missing, 1}, {33, 1, 2}
        };
        VariableType[] types = {VariableType.CONTINUOUS,
            VariableType.BINARY, VariableType.CATEGORICAL};
        MiceOptions options = new MiceOptions(3, 5, 3, 42, 1e-4);
        MiceResult first = MiceImputer.impute(data, types, options);
        MiceResult second = MiceImputer.impute(data, types, options);
        assertEquals(3, first.datasets().size());
        for (int chain = 0; chain < first.datasets().size(); chain++) {
            double[][] completed = first.datasets().get(chain);
            double[][] repeated = second.datasets().get(chain);
            for (int row = 0; row < completed.length; row++) {
                assertArrayEquals(completed[row], repeated[row], 0.0);
                for (int column = 0; column < completed[row].length; column++) {
                    assertFalse(Double.isNaN(completed[row][column]));
                    if (!Double.isNaN(data[row][column]))
                        assertEquals(data[row][column], completed[row][column], 0.0);
                }
            }
            assertTrue(completed[4][2] == 1.0 || completed[4][2] == 2.0);
        }
        assertEquals(1, first.diagnostics().get(0).missingCount());
    }

    @Test void rubinPoolingMatchesHandCalculationAndFiniteDfCorrection() {
        RubinPooling.Estimate pooled = RubinPooling.pool(
            new double[] {1, 2, 3}, new double[] {4, 4, 4}, 20);
        assertEquals(2.0, pooled.estimate(), 0.0);
        assertEquals(4.0, pooled.withinVariance(), 0.0);
        assertEquals(1.0, pooled.betweenVariance(), 0.0);
        assertEquals(16.0 / 3.0, pooled.totalVariance(), 1e-15);
        assertTrue(pooled.degreesOfFreedom() > 0.0
            && pooled.degreesOfFreedom() < 20.0);
        assertTrue(pooled.fractionMissingInformation() > 0.0
            && pooled.fractionMissingInformation() < 1.0);
    }
}
