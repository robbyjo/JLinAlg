/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class SusieTest {
    @Test
    void summaryIbssFindsTwoIndependentSignals() {
        double[] z = {12.0, 0.2, -0.1, -10.0, 0.3};
        double[][] ld = new double[5][5];
        for (int index = 0; index < 5; index++) ld[index][index] = 1.0;
        SusieOptions options = new SusieOptions(
            2, 200, 1e-8, 0.2, false, 0.95, 0.5);

        SusieResult result = Susie.fitSummary(z, ld, 1_000,
            List.of("a", "b", "c", "d", "e"), options, BackendPolicy.CPU);

        assertTrue(result.converged());
        assertTrue(result.pip()[0] > 0.99);
        assertTrue(result.pip()[3] > 0.99);
        assertEquals(2, result.credibleSets().size());
    }

    @Test
    void individualFitReturnsOriginalScaleCoefficientAndIntercept() {
        double[][] design = new double[80][2];
        double[] response = new double[80];
        for (int row = 0; row < 80; row++) {
            design[row][0] = row / 10.0;
            design[row][1] = (row * 17 % 13) / 13.0;
            response[row] = 3.0 + 2.0 * design[row][0];
        }
        SusieResult result = Susie.fit(response, design, List.of("x", "noise"),
            new SusieOptions(1, 200, 1e-8, 10.0, true, 0.95, 0.0),
            BackendPolicy.CPU);
        assertEquals(2.0, result.posteriorMean()[0], 0.05);
        assertEquals(3.0, result.intercept(), 0.2);
        assertTrue(result.pip()[0] > 0.99);
    }
}
