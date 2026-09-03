/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

class PreparedGeeScanTest {
    @Test
    void preparedScanMatchesFullFitsAndParallelExecution() {
        int clusters = 20;
        int size = 4;
        int rows = clusters * size;
        double[] response = new double[rows];
        double[] base = new double[rows * 2];
        double[] predictors = new double[rows * 3];
        int[] id = new int[rows];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < size; visit++) {
                int row = cluster * size + visit;
                base[row * 2] = 1.0;
                base[row * 2 + 1] = visit;
                for (int predictor = 0; predictor < 3; predictor++) {
                    predictors[row * 3 + predictor] =
                        Math.sin((row + 2 * predictor + 1) * 0.37);
                }
                response[row] = 1.2 + 0.3 * visit
                    + 0.2 * predictors[row * 3]
                    + (cluster % 5 - 2) * 0.07;
                id[row] = cluster;
            }
        }
        GeeOptions options = GeeOptions.builder()
            .correlation(GeeCorrelation.EXCHANGEABLE).build();
        try (PreparedGeeScan scan = PreparedGeeScan.prepare(response, base,
                rows, 2, id, null, null, null, GlmFamilies.gaussian(),
                options, BackendPolicy.CPU)) {
            GeeScanResult sequential = scan.scan(predictors, 3,
                List.of("a", "b", "c"), 1);
            GeeScanResult parallel = scan.scan(predictors, 3,
                List.of("a", "b", "c"), 3);
            assertArrayEquals(sequential.coefficients(),
                parallel.coefficients(), 1e-10);
            assertArrayEquals(sequential.standardErrors(),
                parallel.standardErrors(), 1e-10);
            assertTrue(java.util.Arrays.stream(sequential.iterations())
                .allMatch(value -> value > 0));

            double[] fullDesign = new double[rows * 3];
            for (int row = 0; row < rows; row++) {
                fullDesign[row * 3] = base[row * 2];
                fullDesign[row * 3 + 1] = base[row * 2 + 1];
                fullDesign[row * 3 + 2] = predictors[row * 3];
            }
            GeeResult full = Gee.fit(response, fullDesign, rows, 3,
                id, null, GlmFamilies.gaussian(), null, null,
                options, BackendPolicy.CPU);
            assertEquals(full.coefficients()[2],
                sequential.coefficients()[0], 1e-9);
            assertEquals(full.standardErrors()[2],
                sequential.standardErrors()[0], 1e-9);
        }
    }
}
