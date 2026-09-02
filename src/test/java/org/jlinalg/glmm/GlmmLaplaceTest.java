/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

final class GlmmLaplaceTest {
    @Test
    void binomialRandomInterceptUsesLaplaceMarginalLikelihood() {
        int groups = 20;
        int perGroup = 15;
        int observations = groups * perGroup;
        double[] response = new double[observations];
        double[] fixed = new double[observations * 2];
        List<Integer> group = new ArrayList<>(observations);
        long state = 1_234_567L;
        for (int row = 0; row < observations; row++) {
            int cluster = row / perGroup;
            double x = -1.0 + 2.0 * (row % perGroup) / (perGroup - 1.0);
            double random = 1.1 * Math.sin(1.3 * cluster);
            double eta = -0.2 + 0.8 * x + random;
            double probability = 1.0 / (1.0 + Math.exp(-eta));
            state = 48_271L * state % 2_147_483_647L;
            double uniform = state / 2_147_483_647.0;
            response[row] = uniform < probability ? 1.0 : 0.0;
            fixed[row * 2] = 1.0;
            fixed[row * 2 + 1] = x;
            group.add(cluster);
        }
        GlmmLaplaceOptions controls = new GlmmLaplaceOptions(
            14, 60, 1e-4, 1.0, 1e-6, 100.0, null);
        GlmmLaplaceResult fit = GlmmLaplace.fit(response, fixed,
            observations, 2, GlmFamilies.binomial(),
            List.of(VarianceComponent.randomIntercept("group", group)),
            null, null, controls, BackendPolicy.CPU);
        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
        assertEquals(0.8, fit.beta()[1], 0.45);
        assertTrue(fit.varianceComponents()[0] > 1e-6);
        assertEquals(observations, fit.componentPredictor("group").length);
    }
}
