/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class PreparedGamTest {
    @Test
    void warmRefitReusesCompiledSmoothStructure() {
        int observations = 64;
        double[] x = new double[observations];
        double[] first = new double[observations];
        double[] second = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            intercept[row][0] = 1.0;
            first[row] = 2.0 + Math.sin(2.0 * Math.PI * x[row])
                + 0.04 * Math.cos(13.0 * row);
            second[row] = -1.0 + 0.6 * Math.cos(2.0 * Math.PI * x[row])
                + 0.04 * Math.sin(11.0 * row);
        }
        PreparedGam prepared = new PreparedGam(intercept,
            List.of(PSplineTerm.of("s(x)", x, 9)),
            RemlOptions.defaults(), BackendPolicy.CPU);

        GamResult firstFit = prepared.fit(first);
        GamResult secondFit = prepared.refit(firstFit, second);

        assertTrue(firstFit.mixedModel().reml().converged());
        assertTrue(secondFit.mixedModel().reml().converged());
        assertEquals(2.0, firstFit.parametricCoefficients()[0], 0.02);
        assertEquals(-1.0, secondFit.parametricCoefficients()[0], 0.02);
    }
}
