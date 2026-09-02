/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

final class ArmaCorrelationTest {
    @Test
    void ar1MatchesClosedForm() {
        double[] correlation = ArmaCorrelation.ar1(8, 0.65);
        for (int row = 0; row < 8; row++) {
            for (int column = 0; column < 8; column++) {
                assertEquals(Math.pow(0.65, Math.abs(row - column)),
                    correlation[row * 8 + column], 1e-10);
            }
        }
    }
}
