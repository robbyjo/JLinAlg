/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.GaussianSmoothSelectionResult;
import org.jlinalg.gam.SmoothingSelectionOptions;
import org.junit.jupiter.api.Test;

final class AdvancedGamFormulaTest {
    @Test
    void compilesAndFitsTeFormulaOnce() {
        int observations = 100;
        double[] x = new double[observations];
        double[] z = new double[observations];
        double[] y = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = (row % 10) / 9.0;
            z[row] = (row / 10) / 9.0;
            y[row] = 0.5 + Math.sin(2.0 * Math.PI * x[row])
                + 0.3 * Math.cos(2.0 * Math.PI * z[row]);
        }
        ModelTable table = ModelTable.builder(observations)
            .numeric("y", y).numeric("x", x).numeric("z", z).build();
        CompiledQuadraticGamFormula compiled = AdvancedGamFormula.compile(
            "y ~ te(x, z, kx=6, kz=5)", table);
        assertEquals(1, compiled.smoothTerms().size());
        assertEquals(2, compiled.smoothTerms().get(0).penaltyCount());
        GaussianSmoothSelectionResult fit = compiled.fitGaussian(
            new SmoothingSelectionOptions(org.jlinalg.gam.SmoothingCriterion.GCV,
                5, 2.0, -8.0, 10.0, 1e-2, Double.NaN),
            BackendPolicy.CPU);
        assertTrue(fit.effectiveDegreesOfFreedom() > 3.0);
    }

    @Test
    void tiBasisIsConstrainedAgainstMarginalLinearTrends() {
        int observations = 64;
        double[] x = new double[observations];
        double[] z = new double[observations];
        double[] y = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = (row % 8) / 7.0;
            z[row] = (row / 8) / 7.0;
            y[row] = x[row] * z[row];
        }
        ModelTable table = ModelTable.builder(observations)
            .numeric("y", y).numeric("x", x).numeric("z", z).build();
        CompiledQuadraticGamFormula compiled = AdvancedGamFormula.compile(
            "y ~ x + z + ti(x,z,k=5)", table);
        double[] design = compiled.smoothTerms().get(0).design();
        int columns = compiled.smoothTerms().get(0).columns();
        for (int column = 0; column < columns; column++) {
            double sum = 0.0;
            double withX = 0.0;
            double withZ = 0.0;
            for (int row = 0; row < observations; row++) {
                double value = design[row * columns + column];
                sum += value;
                withX += value * x[row];
                withZ += value * z[row];
            }
            assertEquals(0.0, sum, 1e-10);
            assertEquals(0.0, withX, 1e-10);
            assertEquals(0.0, withZ, 1e-10);
        }
    }
}
