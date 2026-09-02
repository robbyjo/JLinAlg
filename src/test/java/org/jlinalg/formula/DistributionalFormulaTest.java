/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.DistributionalFamilies;
import org.jlinalg.distributional.DistributionalOptions;
import org.jlinalg.distributional.DistributionalResult;
import org.junit.jupiter.api.Test;

final class DistributionalFormulaTest {
    @Test
    void meanAndScaleUseSeparateRLikeFormulas() {
        int observations = 120;
        double[] x = new double[observations];
        double[] y = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            double scale = Math.exp(-1.0 + 0.5 * x[row]);
            y[row] = Math.sin(2.0 * Math.PI * x[row])
                + scale * Math.sqrt(2.0) * Math.sin(23.0 * row);
        }
        ModelTable table = ModelTable.builder(observations)
            .numeric("y", y).numeric("x", x).build();
        DistributionalResult fit = DistributionalFormula.fit(
            List.of("y ~ s(x,k=8)", "y ~ x"),
            List.of(List.of(new double[] {0.8}), List.of()), table,
            DistributionalFamilies.gaussianLocationScale(),
            FormulaOptions.defaults(), DistributionalOptions.defaults(),
            BackendPolicy.CPU);
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertTrue(fit.parameter("mu").effectiveDegreesOfFreedom() > 3.0);
        assertTrue(fit.parameter("sigma").coefficients()[1] > 0.2);
    }
}
