/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class MixedFormulaTest {
    @Test
    void compilesSparseRandomInterceptWithoutRuntimeParsing() {
        ModelTable table = ModelTable.builder(9)
            .numeric("y", 0, 1, 2, 4, 5, 6, 8, 9, 10)
            .categorical("subject", "a", "a", "a", "b", "b", "b",
                "c", "c", "c")
            .build();
        CompiledMixedFormula model = MixedFormula.compile(
            "y ~ 1 + (1 | subject)", table);

        assertTrue(model.randomEffects().get(0).sparse());
        assertEquals(3, model.randomEffects().get(0).coefficients());
        assertEquals(5.0, model.fit(
            RemlOptions.builder().initialVariances(10, 2).build(),
            BackendPolicy.CPU).reml().fixedEffects()[0], 1e-10);
    }

    @Test
    void expandsIndependentDoubleBarAndNestedGroupingShorthand() {
        ModelTable table = ModelTable.builder(8)
            .numeric("y", 1, 2, 3, 4, 5, 6, 7, 8)
            .numeric("x", 0, 1, 0, 1, 0, 1, 0, 1)
            .categorical("site", "a", "a", "a", "a", "b", "b", "b", "b")
            .categorical("subject", "u", "u", "v", "v", "u", "u", "v", "v")
            .build();

        CompiledMixedFormula independent = MixedFormula.compile(
            "y ~ x + (1 + x || site)", table);
        assertEquals(2, independent.randomEffects().size());
        assertEquals("1|site", independent.randomEffects().get(0).name());
        assertEquals("0+x|site", independent.randomEffects().get(1).name());

        CompiledMixedFormula nested = MixedFormula.compile(
            "y ~ x + (1 | site/subject)", table);
        assertEquals(2, nested.randomEffects().size());
        assertEquals("1|site:subject", nested.randomEffects().get(1).name());
        assertTrue(nested.fitSparse(
            RemlOptions.builder().initialVariances(1, 1, 1).build(),
            BackendPolicy.CPU).randomCoefficientCount() > 0);
    }

    @Test
    void compilesAndFitsCorrelatedInterceptSlopeBlock() {
        ModelTable table = ModelTable.builder(12)
            .numeric("y", 1.0, 2.0, 2.8, 4.1, 2.1, 3.4,
                4.2, 5.6, 0.4, 1.8, 2.4, 3.7)
            .numeric("x", 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3)
            .categorical("subject", "a", "a", "a", "a",
                "b", "b", "b", "b", "c", "c", "c", "c")
            .build();
        CompiledMixedFormula model = MixedFormula.compile(
            "y ~ x + (1 + x | subject)", table);

        assertEquals(1, model.correlatedRandomEffects().size());
        var fit = model.fitCorrelated(
            RemlOptions.builder().maximumIterations(100).build(),
            BackendPolicy.CPU);
        assertEquals(2, fit.beta().length);
        assertEquals(4, fit.randomEffects().get(0).covariance().length);
        assertTrue(Double.isFinite(fit.randomEffects().get(0)
            .correlation(0, 1)));
    }
}
