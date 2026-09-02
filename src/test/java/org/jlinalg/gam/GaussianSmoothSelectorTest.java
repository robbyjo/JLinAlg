/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

final class GaussianSmoothSelectorTest {
    @Test
    void gcvSelectsFiniteTensorSmoothingAndImprovesOnIntercept() {
        int observations = 100;
        double[] x = new double[observations];
        double[] z = new double[observations];
        double[] y = new double[observations];
        double mean = 0.0;
        for (int row = 0; row < observations; row++) {
            x[row] = (row % 10) / 9.0;
            z[row] = (row / 10) / 9.0;
            y[row] = 0.4 + Math.sin(2.0 * Math.PI * x[row])
                + 0.3 * Math.cos(2.0 * Math.PI * z[row])
                + 0.08 * Math.sin(19.0 * row);
            mean += y[row];
        }
        mean /= observations;
        QuadraticSmoothTerm surface = TensorProductPSplineTerm.of(
            "te(x,z)", x, z, 6, 6);
        SmoothingSelectionOptions controls = new SmoothingSelectionOptions(
            SmoothingCriterion.GCV, 8, 2.0, -10.0, 12.0, 1e-3, Double.NaN);
        GaussianSmoothSelectionResult fit = GaussianSmoothSelector.fit(
            y, intercept(observations), List.of(surface), null,
            controls, BackendPolicy.CPU);
        assertEquals(2, fit.smoothingParameters().get(0).length);
        assertTrue(fit.smoothingParameters().get(0)[0] > 0.0);
        assertTrue(fit.smoothingParameters().get(0)[1] > 0.0);
        assertTrue(fit.effectiveDegreesOfFreedom() > 3.0);
        double nullRss = 0.0;
        for (double value : y) nullRss += (value - mean) * (value - mean);
        double rss = 0.0;
        for (double value : fit.residuals()) rss += value * value;
        assertTrue(rss < 0.2 * nullRss);
    }

    @Test
    void fixedSmoothingDoesNotRunOuterSearch() {
        int observations = 60;
        double[] x = new double[observations];
        double[] y = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = row / 59.0;
            y[row] = 2.0 + Math.sin(2.0 * Math.PI * x[row]);
        }
        GaussianSmoothSelectionResult fit = GaussianSmoothSelector.fitFixed(
            y, intercept(observations),
            List.of(QuadraticSmoothTerm.from(PSplineTerm.of("s(x)", x, 9))),
            List.of(new double[] {0.7}), BackendPolicy.CPU);
        assertEquals(1, fit.evaluations());
        assertEquals(0.7, fit.smoothingParameters().get(0)[0], 0.0);
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
