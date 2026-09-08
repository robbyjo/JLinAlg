/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

final class SemExtensionsTest {
    @Test
    void fimlLatentMeansAndInferenceAreAvailable() {
        double[][] data = new double[60][3];
        for (int row = 0; row < data.length; row++) { double factor = Math.sin(row * 0.2); data[row][0] = factor + 0.1 * Math.cos(row); data[row][1] = 0.7 * factor + 0.1 * Math.sin(row); data[row][2] = -0.4 * factor + 0.1 * Math.cos(row * 0.4); if (row % 7 == 0) data[row][1] = Double.NaN; }
        SemModel model = SemModel.builder("x", "y", "z").regression("y~x", "y", "x", 0.1).build();
        SemFimlResult fiml = SemFiml.fit(data, model);
        LatentMeasurement.Result latent = LatentMeasurement.fit(new double[][] {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {4, 5, 6}, {5, 6, 7}, {6, 7, 8}}, 1);
        SemMeanStructure.Result means = SemMeanStructure.fit(new double[][] {{1, 2}, {2, 3}, {3, 4}});
        SemInference.IndirectEffect indirect = SemInference.indirect(0.5, 0.4, 0.01, 0.02, 0.0, 0.95);
        assertTrue(fiml.converged() || fiml.iterations() > 0); assertEquals(1, latent.factorCount()); assertEquals(2, means.intercepts().length); assertEquals(0.2, indirect.estimate(), 1e-12);
        assertEquals(2, SemInference.ordinalThresholds(new int[] {0, 0, 1, 1, 2, 2}, 3).length);
    }
}
