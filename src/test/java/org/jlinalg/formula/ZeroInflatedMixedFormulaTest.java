/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.ZeroInflatedMixedOptions;
import org.jlinalg.distributional.ZeroInflatedMixedResult;
import org.junit.jupiter.api.Test;

final class ZeroInflatedMixedFormulaTest {
    @Test
    void compilesSeparateMixedPredictorsOnce() {
        int rows = 80;
        double[] y = new double[rows];
        double[] x = new double[rows];
        String[] group = new String[rows];
        for (int row = 0; row < rows; row++) {
            x[row] = (row % 10 - 4.5) / 4.5;
            group[row] = "g" + row / 10;
            y[row] = row % (row / 10 % 2 == 0 ? 4 : 9) == 0
                ? 0.0 : 1.0 + row % 4;
        }
        ModelTable table = ModelTable.builder(rows)
            .numeric("y", y).numeric("x", x)
            .categorical("group", group).build();
        CompiledZeroInflatedMixedFormula formula =
            ZeroInflatedMixedFormula.compilePoisson(
                "y ~ x + (1 || group)",
                "y ~ 1 + (1 || group)", table);
        assertEquals(2, formula.count().columns());
        assertEquals(1, formula.zero().columns());
        assertEquals(1, formula.countRandomEffects().size());
        assertEquals(1, formula.zeroRandomEffects().size());
        ZeroInflatedMixedResult fit = formula.fitPoisson(
            new ZeroInflatedMixedOptions(200, 80, 1e-6, 0.4,
                1e-6, 100.0, 1e-4, 1e4, 20.0, null),
            BackendPolicy.CPU);
        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
    }
}
