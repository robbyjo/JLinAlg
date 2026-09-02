/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

class GeeStressTest {
    @Test
    void singletonAndLargeUnequalClustersFit() {
        List<Double> responses = new ArrayList<>();
        List<double[]> designs = new ArrayList<>();
        List<Integer> ids = new ArrayList<>();
        List<Integer> waves = new ArrayList<>();
        for (int cluster = 0; cluster < 18; cluster++) {
            int size = cluster < 5 ? 1 : cluster == 5 ? 64 : 4;
            for (int visit = 0; visit < size; visit++) {
                responses.add(0.7 + 0.03 * visit
                    + (cluster % 4 - 1.5) * 0.06
                    + Math.sin((cluster + 1.0) * (visit + 1.0)) * 0.01);
                designs.add(new double[] {1.0, visit});
                ids.add(cluster);
                waves.add(visit);
            }
        }
        int n = responses.size();
        double[] response = new double[n];
        double[][] design = new double[n][];
        int[] id = new int[n];
        int[] wave = new int[n];
        for (int row = 0; row < n; row++) {
            response[row] = responses.get(row);
            design[row] = designs.get(row);
            id[row] = ids.get(row);
            wave[row] = waves.get(row);
        }
        GeeResult result = Gee.fit(response, design, id, wave,
            GlmFamilies.gaussian(), null, null,
            GeeOptions.builder().correlation(GeeCorrelation.AR1).build(),
            BackendPolicy.CPU);
        assertTrue(Arrays.stream(result.coefficients()).allMatch(Double::isFinite));
        assertTrue(result.minimumClusterSize() == 1);
        assertTrue(result.maximumClusterSize() == 64);
    }

    @Test
    void nearSeparatedBinomialBiasReductionRemainsFinite() {
        int clusters = 36;
        int size = 2;
        double[] response = new double[clusters * size];
        double[][] design = new double[response.length][2];
        int[] id = new int[response.length];
        int[] wave = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < size; visit++) {
                int row = cluster * size + visit;
                double x = cluster < clusters / 2 ? -1.0 : 1.0;
                response[row] = x > 0.0 ? 1.0 : 0.0;
                if (cluster == 2 && visit == 0) response[row] = 1.0;
                if (cluster == clusters - 3 && visit == 1) response[row] = 0.0;
                design[row][0] = 1.0;
                design[row][1] = x;
                id[row] = cluster;
                wave[row] = visit;
            }
        }
        GeeResult result = Gee.fit(response, design, id, wave,
            GlmFamilies.binomial(), null, null,
            GeeOptions.builder().method(GeeMethod.BIAS_REDUCED)
                .correlation(GeeCorrelation.EXCHANGEABLE).build(),
            BackendPolicy.CPU);
        assertTrue(Arrays.stream(result.coefficients()).allMatch(Double::isFinite));
        assertTrue(Double.isFinite(result.convergenceDiagnostics().estimatingScoreNorm()));
    }

    @Test
    void rankDeficiencyIsRejected() {
        double[] response = {1, 2, 2, 3, 3, 4, 4, 5};
        double[][] design = new double[response.length][3];
        int[] id = {0, 0, 1, 1, 2, 2, 3, 3};
        int[] wave = {0, 1, 0, 1, 0, 1, 0, 1};
        for (int row = 0; row < response.length; row++) {
            design[row][0] = 1.0;
            design[row][1] = row;
            design[row][2] = row;
        }
        assertThrows(IllegalArgumentException.class, () -> Gee.fit(
            response, design, id, wave, GlmFamilies.gaussian(), null, null,
            GeeOptions.defaults(), BackendPolicy.CPU));
    }
}
