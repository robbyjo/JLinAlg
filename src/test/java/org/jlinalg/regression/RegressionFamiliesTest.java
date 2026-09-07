/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

final class RegressionFamiliesTest {
    @Test
    void multivariateOlsRecoversSharedDesign() {
        double[][] x = {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
        double[][] y = new double[6][2];
        for (int row = 0; row < y.length; row++) { y[row][0] = 2 + 3 * row; y[row][1] = -1 + 0.5 * row; }
        MultivariateRegressionResult result = MultivariateRegression.fit(y, x);
        assertEquals(2.0, result.coefficients()[0], 1e-10);
        assertEquals(3.0, result.coefficients()[1], 1e-10);
        assertEquals(-1.0, result.coefficients()[2], 1e-10);
        assertEquals(0.5, result.coefficients()[3], 1e-10);
    }

    @Test
    void multinomialRegressionProducesNormalizedProbabilities() {
        double[][] x = {{1, -2}, {1, -1}, {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
        int[] y = {0, 0, 0, 1, 1, 2, 2, 2};
        MultinomialRegressionResult result = MultinomialRegression.fit(y, x, 3);
        for (int row = 0; row < y.length; row++) {
            double total = 0.0;
            for (int cls = 0; cls < 3; cls++) total += result.probabilities()[row * 3 + cls];
            assertEquals(1.0, total, 1e-12);
        }
        assertTrue(result.logLikelihood() < 0.0);
    }

    @Test
    void quantileKernelSuperAndSemiparametricFitsAreFinite() {
        double[] z = new double[32], y = new double[32];
        double[][] x = new double[32][2];
        for (int row = 0; row < z.length; row++) {
            z[row] = row / 8.0; x[row][0] = 1.0; x[row][1] = z[row];
            y[row] = 1.0 + 0.5 * z[row] + Math.sin(z[row]);
        }
        QuantileRegressionResult quantile = QuantileRegression.fit(y, x, 0.5);
        KernelRegression.Result kernel = KernelRegression.fit(z, y);
        SuperSmoother.Result smoother = SuperSmoother.fit(z, y);
        PartiallyLinearRegression.Result semiparametric = PartiallyLinearRegression.fit(y, x, z, 0.5);
        assertFinite(quantile.coefficients()); assertFinite(kernel.fittedValues());
        assertFinite(smoother.fittedValues()); assertFinite(semiparametric.fittedValues());
        assertTrue(quantile.objective() >= 0.0);
    }

    private static void assertFinite(double[] values) {
        for (double value : values) assertTrue(Double.isFinite(value));
    }
}
