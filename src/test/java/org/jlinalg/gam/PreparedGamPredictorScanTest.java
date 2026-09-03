/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

final class PreparedGamPredictorScanTest {
    @Test
    void retainedContextMatchesOrdinaryGaussianGam() {
        int observations = 70;
        double[] response = new double[observations];
        double[] x = new double[observations];
        double[] fixed = new double[observations];
        java.util.Arrays.fill(fixed, 1.0);
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            response[row] = 1.2 + Math.sin(2.0 * Math.PI * x[row])
                + 0.03 * Math.cos(13.0 * row);
        }
        GamResult ordinary = Gam.fitGaussian(response, fixed, observations, 1,
            java.util.List.of(PSplineTerm.of("s(x)", x, 10)),
            RemlOptions.defaults(), BackendPolicy.CPU);
        PreparedGamPredictorScan prepared = new PreparedGamPredictorScan(
            response, fixed, observations, 1,
            RemlOptions.defaults(), BackendPolicy.CPU);
        GamResult retained = prepared.fit("s(x)", x, 10);

        assertArrayEquals(ordinary.fittedValues(), retained.fittedValues(), 1e-10);
        assertEquals(ordinary.smoothTerms().get(0).effectiveDegreesOfFreedom(),
            retained.smoothTerms().get(0).effectiveDegreesOfFreedom(), 1e-10);
        prepared.close();
        assertThrows(IllegalStateException.class,
            () -> prepared.fit("s(x)", x, 10));
    }
}
