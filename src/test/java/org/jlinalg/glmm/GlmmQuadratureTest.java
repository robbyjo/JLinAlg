/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

final class GlmmQuadratureTest {
    @Test
    void quadratureFitsRandomInterceptBinaryData() {
        double[] y = {0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1}; double[][] x = new double[y.length][2]; List<String> groups = List.of("a", "a", "a", "b", "b", "b", "c", "c", "c", "d", "d", "d");
        for (int row = 0; row < y.length; row++) { x[row][0] = 1.0; x[row][1] = row % 3; }
        GlmmQuadratureResult result = GlmmQuadrature.fit(y, x, groups, GlmFamilies.binomial());
        assertEquals(10, result.nodes()); assertEquals(y.length, result.fittedMeans().length); assertTrue(Double.isFinite(result.logLikelihood()));
    }
}
