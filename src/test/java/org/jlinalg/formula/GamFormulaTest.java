/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GamResult;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class GamFormulaTest {
    @Test
    void compilesFixedTermsOffsetsAndSmoothOptionsOnce() {
        int observations = 60;
        double[] y = new double[observations];
        double[] x = new double[observations];
        double[] exposure = new double[observations];
        String[] group = new String[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            exposure[row] = 0.25;
            group[row] = (row & 1) == 0 ? "A" : "B";
            y[row] = 2.0 + exposure[row]
                + ("B".equals(group[row]) ? 0.5 : 0.0)
                + Math.sin(2.0 * Math.PI * x[row])
                + 0.03 * Math.cos(11.0 * row);
        }
        ModelTable table = ModelTable.builder(observations)
            .numeric("y", y)
            .numeric("x", x)
            .numeric("exposure", exposure)
            .categorical("group", group)
            .build();

        CompiledGamFormula formula = GamFormula.compile(
            "y ~ group + s(x, bs='ps', k=9, degree=3, m=2)"
                + " + offset(exposure)", table);
        GamResult result = formula.fitGaussian(
            RemlOptions.defaults(), BackendPolicy.CPU);

        assertEquals(2, formula.parametricColumns());
        assertEquals(1, formula.smoothTerms().size());
        assertEquals(9, formula.smoothTerms().get(0).basisDimension());
        assertTrue(result.mixedModel().reml().converged(),
            result.mixedModel().reml().convergenceMessage());
        assertEquals(0.5, result.parametricCoefficients()[1], 0.03);
    }
}
