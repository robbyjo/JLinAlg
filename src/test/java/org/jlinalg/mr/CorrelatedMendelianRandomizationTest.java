/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class CorrelatedMendelianRandomizationTest {
    @Test
    void generalizedIvwWithIdentityMatchesIndependentIvw() {
        List<HarmonizedInstrument> instruments = List.of(
            instrument("a", 0.1, 0.05),
            instrument("b", 0.2, 0.10),
            instrument("c", 0.3, 0.15));
        double[][] identity = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        CorrelatedMrEstimate generalized = CorrelatedMendelianRandomization.ivw(
            instruments, identity, false, 0.95, BackendPolicy.CPU);
        MrEstimate independent = MendelianRandomization.ivw(
            instruments, false, 0.95);

        assertEquals(independent.estimate(),
            generalized.estimate().estimate(), 1e-14);
        assertEquals(independent.standardError(),
            generalized.estimate().standardError(), 1e-14);
        assertEquals("cpu", generalized.backend().selectedBackend());
    }

    @Test
    void generalizedEggerHandlesOrientationAndLdSigns() {
        List<HarmonizedInstrument> instruments = List.of(
            instrument("a", -0.1, -0.07),
            instrument("b", 0.2, 0.12),
            instrument("c", 0.3, 0.17),
            instrument("d", 0.4, 0.22));
        double[][] correlation = equicorrelation(4, 0.2);

        CorrelatedMrEggerResult result = CorrelatedMendelianRandomization.egger(
            instruments, correlation, 0.95, BackendPolicy.CPU);

        assertEquals(0.5, result.estimate().slope().estimate(), 1e-12);
        assertEquals(0.02, result.estimate().intercept(), 1e-13);
    }

    @Test
    void rejectsNonPositiveDefiniteLd() {
        List<HarmonizedInstrument> instruments = List.of(
            instrument("a", 0.1, 0.05),
            instrument("b", 0.2, 0.10),
            instrument("c", 0.3, 0.15));
        assertThrows(IllegalArgumentException.class,
            () -> CorrelatedMendelianRandomization.ivw(
                instruments, equicorrelation(3, 1.0), false,
                0.95, BackendPolicy.CPU));
    }

    private static HarmonizedInstrument instrument(
            String id, double exposure, double outcome) {
        return new HarmonizedInstrument(id, "A", "C",
            exposure, 0.01, outcome, 0.01,
            0.2, 0.2, false, false);
    }

    private static double[][] equicorrelation(int size, double correlation) {
        double[][] result = new double[size][size];
        for (int row = 0; row < size; row++) {
            for (int column = 0; column < size; column++) {
                result[row][column] = row == column ? 1.0 : correlation;
            }
        }
        return result;
    }
}
