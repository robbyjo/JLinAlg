/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class CoxDiagnosticsTest {
    @Test
    void exportsClusterRobustInfluenceAndPhDiagnostics() {
        double[] start = {0, 0, 0, 0, 1.1, 1.2, 1.5, 2.2, 2.3, 3.2};
        double[] stop = {1, 2, 2, 3, 4, 4, 5, 6, 7, 8};
        boolean[] event = {
            true, true, false, true, true, true, false, true, false, true
        };
        double[][] x = {
            {-1, 0}, {-0.5, 1}, {0, 0}, {0.5, 1}, {1, 0},
            {-1, 1}, {-0.5, 0}, {0, 1}, {0.5, 0}, {1, 1}
        };
        CoxSurvivalData survival = new CoxSurvivalData(
            start, stop, event, null);
        CoxResult fit = CoxRegression.fit(survival, x, null,
            CoxOptions.defaults(), BackendPolicy.CPU);
        CoxDiagnosticsResult diagnostics = CoxDiagnostics.analyze(
            survival, x, fit, null,
            List.of("a", "a", "b", "b", "c", "c", "d", "d", "e", "e"));
        assertEquals(5, diagnostics.clusters());
        // survival::coxph(... + cluster(id)), survival 3.8-3 / R 4.6.1.
        assertEquals(-0.6238980, fit.beta()[0], 1e-6);
        assertEquals(0.4842102, fit.beta()[1], 1e-6);
        assertEquals(0.5622623, diagnostics.robustStandardErrors()[0], 1e-6);
        assertEquals(0.6003388, diagnostics.robustStandardErrors()[1], 1e-6);
        assertEquals(10, diagnostics.residuals().martingale().length);
        assertEquals(20, diagnostics.residuals().score().length);
        assertEquals(20, diagnostics.residuals().dfbeta().length);
        assertEquals(7, diagnostics.residuals().schoenfeldRows().length);
        assertEquals(14, diagnostics.residuals().schoenfeld().length);
        assertEquals(2, diagnostics.robustStandardErrors().length);
        for (double value : diagnostics.robustStandardErrors())
            assertTrue(value > 0.0 && Double.isFinite(value));
        assertEquals(7, diagnostics.proportionalHazards().events());
        assertEquals(2,
            diagnostics.proportionalHazards().globalDegreesOfFreedom());
        assertTrue(diagnostics.proportionalHazards().globalPValue() >= 0.0
            && diagnostics.proportionalHazards().globalPValue() <= 1.0);
        double martingaleSum = java.util.Arrays.stream(
            diagnostics.residuals().martingale()).sum();
        assertEquals(0.0, martingaleSum, 1e-10);
        double[] score = diagnostics.residuals().score();
        for (int column = 0; column < 2; column++) {
            double sum = 0.0;
            for (int row = 0; row < 10; row++) sum += score[row * 2 + column];
            assertEquals(0.0, sum, 1e-7);
        }
    }
}
