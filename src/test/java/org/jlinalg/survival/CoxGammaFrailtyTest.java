/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class CoxGammaFrailtyTest {
    @Test
    void fitsSharedGammaFrailtyAtFixedAndProfiledVariance() {
        double[] time = {1, 2, 3, 5, 2, 4, 6, 8, 1.5, 3.5, 7, 9};
        boolean[] event = {true, true, false, true, true, false,
            true, false, true, true, false, true};
        double[][] x = {{-1}, {-0.5}, {0}, {0.5}, {-0.8}, {-0.2},
            {0.3}, {0.9}, {-0.7}, {0.1}, {0.6}, {1.1}};
        List<String> groups = List.of("a", "a", "a", "a",
            "b", "b", "b", "b", "c", "c", "c", "c");
        CoxSurvivalData survival = CoxSurvivalData.rightCensored(time, event);
        CoxGammaFrailtyOptions options = new CoxGammaFrailtyOptions(
            CoxOptions.defaults(), 0.2, 1e-5, 5.0, 20, 1e-3);

        CoxGammaFrailtyResult fixed = CoxGammaFrailty.fitAtVariance(
            survival, x, groups, null, 0.2, options, BackendPolicy.CPU);
        CoxGammaFrailtyResult profiled = CoxGammaFrailty.fit(
            survival, x, groups, null, options, BackendPolicy.CPU);

        assertEquals(0.2, fixed.frailtyVariance(), 0.0);
        // survival::coxph(... frailty(..., method="fixed", theta=.2)), R 4.6.1.
        assertEquals(-5.2414724, fixed.beta()[0], 1e-6);
        assertEquals(2.0247874, fixed.standardErrors()[0], 1e-6);
        assertEquals(0.2078657, fixed.logFrailties()[0], 1e-6);
        assertEquals(-0.4582663, fixed.logFrailties()[1], 1e-6);
        assertEquals(0.1280177, fixed.logFrailties()[2], 1e-6);
        assertEquals(3, fixed.frailties().length);
        assertTrue(java.util.Arrays.stream(fixed.frailties())
            .allMatch(value -> value > 0.0 && Double.isFinite(value)));
        assertTrue(Double.isFinite(fixed.beta()[0]));
        assertTrue(Double.isFinite(fixed.laplaceLogLikelihood()));
        assertTrue(profiled.frailtyVariance() >= options.minimumVariance());
        assertTrue(profiled.frailtyVariance() <= options.maximumVariance());
        assertTrue(Double.isFinite(profiled.laplaceLogLikelihood()));
    }
}
