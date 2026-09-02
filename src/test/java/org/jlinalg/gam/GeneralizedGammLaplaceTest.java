/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.GlmmLaplaceOptions;
import org.junit.jupiter.api.Test;

final class GeneralizedGammLaplaceTest {
    @Test
    void poissonSmoothUsesLaplaceIntegration() {
        int observations = 90;
        double[] x = new double[observations];
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            double mean = Math.exp(0.5 + 0.55 * Math.sin(2.0 * Math.PI * x[row]));
            double uniform = ((row * 41 + 11) % 101) / 101.0;
            response[row] = uniform < Math.exp(-mean) ? 0.0
                : uniform < Math.exp(-mean) * (1.0 + mean) ? 1.0
                : uniform < Math.exp(-mean) * (1.0 + mean + mean * mean / 2.0)
                    ? 2.0 : 3.0;
            intercept[row][0] = 1.0;
        }
        GlmmLaplaceOptions controls = new GlmmLaplaceOptions(
            10, 50, 1e-3, 1.0, 1e-5, 1e4, null);
        GeneralizedGammLaplaceResult fit = GeneralizedGammLaplace.fit(
            response, intercept, List.of(PSplineTerm.of("s(x)", x, 7)),
            GlmFamilies.poisson(), List.of(), null, null,
            controls, BackendPolicy.CPU);
        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
        assertEquals(observations, fit.smoothContributions().get("s(x)").length);
        assertTrue(fit.smoothingParameters().get("s(x)") > 0.0);
        assertTrue(fit.parametricStandardErrors()[0] > 0.0);
    }
}
