/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.junit.jupiter.api.Test;

class GeneralizedGamTest {
    @Test
    void poissonPqlSmoothImprovesOverAnInterceptOnlyGlm() {
        int observations = 120;
        double[] x = sequence(observations, 0.0, 1.0);
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            intercept[row][0] = 1.0;
            double mean = Math.exp(1.2 + 0.85
                * Math.sin(2.0 * Math.PI * x[row]));
            double fraction = ((row * 37) % 101) / 101.0;
            response[row] = Math.floor(mean + fraction);
        }
        double interceptOnlyDeviance = Glm.fit(response, intercept,
            GlmFamilies.poisson(), null, null, GlmOptions.defaults(),
            BackendPolicy.CPU).deviance();

        GeneralizedGamResult result = GeneralizedGam.fit(
            response, intercept, List.of(PSplineTerm.of("s(x)", x, 9)),
            GlmFamilies.poisson(), null, null, GlmmPqlOptions.defaults(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertTrue(result.conditionalDeviance() < 0.35 * interceptOnlyDeviance);
        assertTrue(result.smoothTerms().get(0).effectiveDegreesOfFreedom() > 2.0);
        double[] rowMajorIntercept = new double[observations];
        java.util.Arrays.fill(rowMajorIntercept, 1.0);
        double[] predicted = result.predictMean(
            rowMajorIntercept, observations, List.of(x), null);
        double[] fitted = result.fittedMeans();
        for (int row = 0; row < observations; row++) {
            assertEquals(fitted[row], predicted[row], 2e-4);
        }
    }

    private static double[] sequence(int length, double lower, double upper) {
        double[] result = new double[length];
        for (int index = 0; index < length; index++) {
            result[index] = lower + (upper - lower) * index / (length - 1.0);
        }
        return result;
    }
}
