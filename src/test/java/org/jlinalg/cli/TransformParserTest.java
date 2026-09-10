/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import java.util.List;
import org.jlinalg.pipeline.OmicsTransform;
import org.junit.jupiter.api.Test;

class TransformParserTest {
    private static final double MAD_SCALE = 1.482602218505602;

    @Test
    void parsesDefaultAndParameterizedMadWinsorization() {
        double[] values = {1, 2, 3, 4, 100, Double.NaN};
        OmicsTransform defaults = TransformParser.parse(
            List.of("<omics>=winsor_mad"), List.of());
        OmicsTransform parameterized = TransformParser.parse(
            List.of("<omics> = identity() | winsor_mad(k=2)"), List.of());

        assertArrayEquals(
            new double[] {1, 2, 3, 4, 3 + 4 * MAD_SCALE, Double.NaN},
            defaults.apply(values), 1e-12);
        assertArrayEquals(
            new double[] {1, 2, 3, 4, 3 + 2 * MAD_SCALE, Double.NaN},
            parameterized.apply(values), 1e-12);
    }
}
