/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Random;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class ExactArmaTest {
    @Test
    void exactArOneHandlesMissingValuesAndReturnsCoefficientSe() {
        double[] series = simulateAr(180, 0.65, 2.0, 1234L);
        series[20] = Double.NaN;
        series[77] = Double.NaN;
        ExactArmaResult result = ExactArma.fit(
            series, ArimaOrder.ar(1), true, BackendPolicy.CPU);
        assertEquals(0.65, result.autoregressive()[0], 0.12);
        assertEquals(178, result.observedValues());
        assertEquals(2, result.standardErrors().length);
        assertTrue(result.standardErrors()[0] > 0.0);
    }

    @Test
    void panelLikelihoodSharesDynamicsAcrossIndependentSeries() {
        ExactArmaResult result = ExactArma.fitPanel(List.of(
            simulateAr(90, 0.5, 1.0, 1L),
            simulateAr(100, 0.5, 1.0, 2L)),
            ArimaOrder.ar(1), true, BackendPolicy.CPU);
        assertEquals(2, result.seriesCount());
        assertEquals(0.5, result.autoregressive()[0], 0.15);
    }

    @Test
    void automaticSelectionEvaluatesRequestedGrid() {
        double[] series = simulateAr(120, 0.7, 0.0, 91L);
        ArimaSelectionResult selected = AutomaticArima.select(series, 1, 0, 1);
        assertEquals(4, selected.candidates().size());
        assertTrue(Double.isFinite(selected.bestModel().aicc()));
    }

    private static double[] simulateAr(int length, double ar, double mean, long seed) {
        Random random = new Random(seed);
        double[] all = new double[length + 200];
        for (int index = 0; index < all.length; index++) {
            all[index] = mean + random.nextGaussian()
                + (index == 0 ? 0.0 : ar * (all[index - 1] - mean));
        }
        return java.util.Arrays.copyOfRange(all, 200, all.length);
    }
}
