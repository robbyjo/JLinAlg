/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class CategoricalGeeTest {
    @Test
    void partialProportionalOddsSupportsCorrelationAndLocalOddsRatios() {
        int clusters = 30;
        int visits = 2;
        int[] response = new int[clusters * visits];
        double[][] x = new double[response.length][2];
        int[] id = new int[response.length];
        int[] wave = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < visits; visit++) {
                int row = cluster * visits + visit;
                response[row] = (cluster + visit + cluster / 5) % 3;
                x[row][0] = visit;
                x[row][1] = cluster % 2;
                id[row] = cluster;
                wave[row] = visit;
            }
        }
        GeeOptions options = GeeOptions.builder()
            .correlation(GeeCorrelation.EXCHANGEABLE)
            .association(GeeAssociation.ODDS_RATIO).build();
        OrdinalGeeResult result = OrdinalGee.fitPartial(response, x, id, wave,
            3, new boolean[] {true, false}, options, BackendPolicy.CPU);

        assertEquals(2, result.thresholds().length);
        assertEquals(3, result.coefficients().length);
        assertTrue(Arrays.stream(result.cutoffSpecificCoefficients())
            .allMatch(Double::isFinite));
        assertEquals(GeeAssociation.ODDS_RATIO, result.fit().association());
    }

    @Test
    void nominalMultinomialGeeReturnsNormalizedProbabilities() {
        int clusters = 45;
        int visits = 2;
        int[] response = new int[clusters * visits];
        double[][] design = new double[response.length][2];
        int[] id = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < visits; visit++) {
                int row = cluster * visits + visit;
                response[row] = (cluster + 2 * visit + cluster / 7) % 3;
                design[row][0] = 1.0;
                design[row][1] = visit + (cluster % 5) * 0.1;
                id[row] = cluster;
            }
        }
        NominalGeeResult result = NominalGee.fit(response, design, id, 3,
            GeeOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged());
        assertEquals(4, result.coefficients().length);
        double[] probabilities = result.fittedProbabilities();
        for (int row = 0; row < response.length; row++) {
            double sum = probabilities[row * 3]
                + probabilities[row * 3 + 1] + probabilities[row * 3 + 2];
            assertEquals(1.0, sum, 1e-12);
        }
    }
}
