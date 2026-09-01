/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.inference;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class AssociationStatisticsTest {
    @Test
    void exposesRegularLog10AndNegativeLog10PValues() {
        AssociationStatistics statistics = AssociationStatistics.normal(
            new double[] {2.0}, new double[] {1.0});

        double p = statistics.pValues()[0];
        assertEquals(Math.log10(p), statistics.log10PValues()[0], 1e-14);
        assertEquals(-Math.log10(p),
            statistics.negativeLog10PValues()[0], 1e-14);
        assertArrayEquals(statistics.pValues(),
            statistics.pValues(PValueScale.REGULAR), 0.0);
        assertArrayEquals(statistics.log10PValues(),
            statistics.pValues(PValueScale.LOG10), 0.0);
        assertArrayEquals(statistics.negativeLog10PValues(),
            statistics.pValues(PValueScale.NEGATIVE_LOG10), 0.0);
        assertArrayEquals(statistics.beta(), statistics.effectSizes(), 0.0);
    }

    @Test
    void computesExtremeLogPDirectlyWithoutRegularPUnderflow() {
        AssociationStatistics statistics = AssociationStatistics.normal(
            new double[] {40.0}, new double[] {1.0});

        assertEquals(0.0, statistics.pValues()[0]);
        assertTrue(Double.isFinite(statistics.negativeLog10PValues()[0]));
        assertTrue(statistics.negativeLog10PValues()[0] > 300.0);
    }
}
