/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class GlmmPqlTest {
    @Test
    void poissonRandomInterceptFindsOrderedGroupEffects() {
        double[] response = {
            0.0, 1.0, 1.0,
            1.0, 2.0, 3.0,
            3.0, 4.0, 5.0,
            7.0, 8.0, 9.0
        };
        double[][] fixed = new double[response.length][1];
        for (double[] row : fixed) {
            row[0] = 1.0;
        }
        double[] group = groupedRelationship(response.length, 3);
        GlmmPqlOptions options = GlmmPqlOptions.builder()
            .maximumIterations(50)
            .relativeTolerance(1e-5)
            .remlOptions(RemlOptions.builder()
                .initialVariances(0.5)
                .relativeTolerance(1e-9)
                .scoreTolerance(1e-6)
                .build())
            .build();

        GlmmPqlResult result = GlmmPql.fit(
            response, fixed, GlmFamilies.poisson(),
            List.of(new VarianceComponent("group", response.length, group)),
            null, null, options, BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertTrue(result.varianceComponents()[0] > 0.05);
        assertEquals(response.length - 1 - 1,
            result.degreesOfFreedom()[0]);
        assertTrue(Double.isFinite(result.tStatistics()[0]));
        assertTrue(result.pValues()[0] >= 0.0 && result.pValues()[0] <= 1.0);
        double[] random = result.randomLinearPredictor();
        assertTrue(random[0] < random[3]);
        assertTrue(random[3] < random[6]);
        assertTrue(random[6] < random[9]);
        for (int groupIndex = 0; groupIndex < 4; groupIndex++) {
            int first = groupIndex * 3;
            assertEquals(random[first], random[first + 1], 1e-10);
            assertEquals(random[first], random[first + 2], 1e-10);
        }
    }

    @Test
    void gaussianFamilyDirectsCallerToExactReml() {
        assertThrows(IllegalArgumentException.class,
            () -> GlmmPql.fit(
                new double[] {1.0, 2.0, 3.0},
                new double[][] {{1.0}, {1.0}, {1.0}},
                GlmFamilies.gaussian(),
                List.of(VarianceComponent.identity("random", 3)),
                null, null, GlmmPqlOptions.defaults(), BackendPolicy.CPU));
    }

    @Test
    void binomialRandomInterceptProducesValidProbabilities() {
        double[] response = {
            0, 0, 0, 0, 1,
            0, 0, 0, 1, 0,
            0, 0, 1, 1, 0,
            1, 1, 1, 0, 0,
            1, 1, 1, 1, 0,
            1, 1, 1, 0, 1
        };
        double[][] fixed = new double[response.length][1];
        for (double[] row : fixed) {
            row[0] = 1.0;
        }
        double[] group = groupedRelationship(response.length, 5);
        GlmmPqlOptions options = GlmmPqlOptions.builder()
            .maximumIterations(50)
            .relativeTolerance(1e-5)
            .remlOptions(RemlOptions.builder()
                .initialVariances(0.5)
                .scoreTolerance(1e-6)
                .build())
            .build();

        GlmmPqlResult result = GlmmPql.fit(
            response, fixed, GlmFamilies.binomial(),
            List.of(new VarianceComponent("group", response.length, group)),
            null, null, options, BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertTrue(result.varianceComponents()[0] > 0.01);
        for (double probability : result.fittedMeans()) {
            assertTrue(probability > 0.0 && probability < 1.0);
        }
        assertTrue(result.randomLinearPredictor()[0]
            < result.randomLinearPredictor()[25]);
    }

    private static double[] groupedRelationship(int observations, int groupSize) {
        double[] result = new double[observations * observations];
        for (int row = 0; row < observations; row++) {
            for (int column = 0; column < observations; column++) {
                if (row / groupSize == column / groupSize) {
                    result[row * observations + column] = 1.0;
                }
            }
        }
        return result;
    }
}
