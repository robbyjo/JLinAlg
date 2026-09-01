/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.OlsOptions;
import org.junit.jupiter.api.Test;

class FormulaTest {
    @Test
    void compilesTreatmentContrastsInteractionsOffsetAndWeightsOnce() {
        ModelTable table = ModelTable.builder(6)
            .numeric("y", 11, 13, 15, 17, 19, 21)
            .numeric("x", 0, 1, 2, 3, 4, 5)
            .numeric("exposure", 10, 10, 10, 10, 10, 10)
            .numeric("w", 1, 2, 1, 2, 1, 2)
            .categorical("sex", "F", "M", "F", "M", "F", "M")
            .build();

        CompiledFormula model = Formula.compile(
            "y ~ x * sex + offset(exposure)", table,
            new FormulaOptions(ContrastCoding.TREATMENT, "w"));

        assertEquals(List.of("(Intercept)", "x", "sexM", "x:sexM"),
            model.coefficientNames());
        assertArrayEquals(new double[] {
            1, 0, 0, 0,
            1, 1, 1, 1,
            1, 2, 0, 0,
            1, 3, 1, 3,
            1, 4, 0, 0,
            1, 5, 1, 5
        }, model.design());
        assertEquals(2.0, model.fitOls(OlsOptions.defaults(), BackendPolicy.CPU)
            .coefficients()[1], 1e-12);
    }
}
