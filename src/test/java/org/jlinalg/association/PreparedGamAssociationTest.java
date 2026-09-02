/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GaussianSmoothSelectionResult;
import org.jlinalg.gam.GaussianSmoothSelector;
import org.jlinalg.gam.PSplineTerm;
import org.jlinalg.gam.QuadraticSmoothTerm;
import org.junit.jupiter.api.Test;

final class PreparedGamAssociationTest {
    @Test
    void retainedNullProjectionFindsMarkerEffectInParallel() {
        int observations = 180;
        double[] x = new double[observations];
        double[] response = new double[observations];
        double[][] markers = new double[observations][2];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            markers[row][0] = Math.sqrt(2.0) * Math.sin(31.0 * row + 0.2);
            markers[row][1] = Math.sqrt(2.0) * Math.cos(43.0 * row + 0.7);
            response[row] = 1.0 + Math.sin(2.0 * Math.PI * x[row])
                + 0.7 * markers[row][0] + 0.08 * Math.sin(17.0 * row);
        }
        GaussianSmoothSelectionResult nullFit = GaussianSmoothSelector.fitFixed(
            response, intercept(observations),
            List.of(QuadraticSmoothTerm.from(
                PSplineTerm.of("s(x)", x, 10))),
            List.of(new double[] {0.5}), BackendPolicy.CPU);
        AssociationEngineOptions options = AssociationEngineOptions.cpuParallel()
            .withParallelism(2);
        try (PreparedGamAssociation scan =
                new PreparedGamAssociation(nullFit, options)) {
            AssociationBatchResult result = scan.scan(
                markers, List.of("causal", "null"));
            assertTrue(result.successful());
            assertEquals(0.7, result.beta()[0], 0.03);
            assertTrue(Math.abs(result.statistics()[0]) > 10.0);
            assertTrue(Math.abs(result.beta()[1]) < 0.04);
            assertTrue(result.negativeLog10PValues()[0] > 10.0);
        }
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
