/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.Test;

class MendelianRandomizationTest {
    @Test
    void ivwAndWaldRecoverCommonRatio() {
        List<HarmonizedInstrument> instruments = commonRatioInstruments();

        MrEstimate fixed = MendelianRandomization.ivw(
            instruments, false, 0.95);
        MrEstimate random = MendelianRandomization.ivw(
            instruments, true, 0.95);
        List<WaldRatio> ratios = MendelianRandomization.waldRatios(instruments);

        assertEquals(0.5, fixed.estimate(), 1e-14);
        assertEquals(Math.sqrt(1.0 / 3000.0), fixed.standardError(), 1e-14);
        assertEquals(0.0, fixed.cochranQ(), 1e-25);
        assertEquals(fixed.standardError(), random.standardError(), 1e-14);
        assertEquals(0.5, ratios.get(0).estimate(), 1e-14);
        assertEquals(100.0, ratios.get(0).fStatistic(), 1e-12);
    }

    @Test
    void eggerRecoversSlopeAndDirectionalIntercept() {
        List<HarmonizedInstrument> instruments = List.of(
            instrument("a", 0.1, 0.07),
            instrument("b", 0.2, 0.12),
            instrument("c", 0.3, 0.17),
            instrument("d", 0.4, 0.22));

        MrEggerResult result = MendelianRandomization.egger(instruments, 0.95);

        assertEquals(0.5, result.slope().estimate(), 1e-13);
        assertEquals(0.02, result.intercept(), 1e-14);
        assertEquals(Math.sqrt(0.002),
            result.slope().standardError(), 1e-13);
        assertTrue(result.iSquaredGx() > 0.99);
    }

    @Test
    void fullAnalysisIncludesMedianStrengthAndLeaveOneOut() {
        MrAnalysisResult result = MendelianRandomization.analyze(
            commonRatioInstruments(), new MrOptions(0.95, 200, 42L));

        assertEquals(0.5, result.weightedMedian().estimate(), 1e-14);
        assertTrue(result.weightedMedian().standardError() > 0.0);
        assertEquals(750.0, result.meanFStatistic(), 1e-10);
        assertEquals(4, result.leaveOneOut().size());
        for (LeaveOneOutEstimate estimate : result.leaveOneOut()) {
            assertEquals(0.5, estimate.estimate().estimate(), 1e-14);
        }
    }

    @Test
    void rejectsDuplicateAndInsufficientInstruments() {
        HarmonizedInstrument instrument = instrument("a", 0.1, 0.05);
        assertThrows(IllegalArgumentException.class,
            () -> MendelianRandomization.analyze(List.of(instrument, instrument, instrument)));
        assertThrows(IllegalArgumentException.class,
            () -> MendelianRandomization.ivw(List.of(instrument), false, 0.95));
    }

    private static List<HarmonizedInstrument> commonRatioInstruments() {
        return List.of(
            instrument("a", 0.1, 0.05),
            instrument("b", 0.2, 0.10),
            instrument("c", 0.3, 0.15),
            instrument("d", 0.4, 0.20));
    }

    private static HarmonizedInstrument instrument(
            String id, double exposureEffect, double outcomeEffect) {
        return new HarmonizedInstrument(id, "A", "C",
            exposureEffect, 0.01, outcomeEffect, 0.01,
            0.2, 0.2, false, false);
    }
}
