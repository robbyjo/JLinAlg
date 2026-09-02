/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

/** Cross-language regression test against lme4's nAGQ=1 Laplace fit. */
final class Lme4LaplaceReferenceTest {
    @Test
    void matchesLme4BinomialRandomInterceptReference() throws IOException {
        Properties reference = new Properties();
        try (InputStream input = getClass().getResourceAsStream(
                "/r-reference/lme4-binomial-laplace.properties")) {
            if (input == null) throw new IOException("missing lme4 reference fixture");
            reference.load(input);
        }

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
            double eta = -0.2 + 0.8 * x + 1.1 * Math.sin(1.3 * cluster);
            double probability = 1.0 / (1.0 + Math.exp(-eta));
            state = 48_271L * state % 2_147_483_647L;
            response[row] = state / 2_147_483_647.0 < probability ? 1.0 : 0.0;
            fixed[row * 2] = 1.0;
            fixed[row * 2 + 1] = x;
            group.add(cluster);
        }

        GlmmLaplaceResult fit = GlmmLaplace.fit(response, fixed,
            observations, 2, GlmFamilies.binomial(),
            List.of(VarianceComponent.randomIntercept("group", group)),
            null, null,
            new GlmmLaplaceOptions(30, 100, 1e-7, 1.0, 1e-8, 100.0, null),
            BackendPolicy.CPU);

        assertTrue(fit.converged());
        assertEquals(value(reference, "intercept"), fit.beta()[0], 0.01);
        assertEquals(value(reference, "slope"), fit.beta()[1], 0.04);
        assertEquals(value(reference, "randomVariance"),
            fit.varianceComponents()[0], 0.02);
        assertEquals(value(reference, "logLikelihood"),
            fit.marginalLogLikelihood(), 0.05);
        assertEquals(value(reference, "fitted0"), fit.fittedMeans()[0], 0.02);
        assertEquals(value(reference, "fitted149"), fit.fittedMeans()[149], 0.02);
        assertEquals(value(reference, "fitted299"), fit.fittedMeans()[299], 0.02);
    }

    private static double value(Properties properties, String key) {
        return Double.parseDouble(properties.getProperty(key));
    }
}
