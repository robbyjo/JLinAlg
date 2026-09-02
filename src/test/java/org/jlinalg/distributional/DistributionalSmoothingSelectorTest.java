/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PSplineTerm;
import org.jlinalg.gam.QuadraticSmoothTerm;
import org.junit.jupiter.api.Test;

final class DistributionalSmoothingSelectorTest {
    @Test
    void selectsMeanSmoothWhileScaleStaysParametric() {
        int observations = 90;
        double[] x = new double[observations];
        double[] y = new double[observations];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            y[row] = Math.sin(2.0 * Math.PI * x[row])
                + 0.25 * Math.sin(29.0 * row);
        }
        double[][] intercept = intercept(observations);
        DistributionalSmoothingResult result = DistributionalSmoothingSelector.fit(
            y, List.of(intercept, intercept),
            List.of(List.of(QuadraticSmoothTerm.from(
                    PSplineTerm.of("s(x)", x, 8))), List.of()),
            null, DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(),
            new DistributionalSmoothingOptions(5, 2.0, -8.0, 10.0, 1e-2),
            BackendPolicy.CPU);
        assertTrue(result.fit().converged(), result.fit().convergenceMessage());
        assertTrue(result.smoothingParameters().get(0).get(0)[0] > 0.0);
        assertTrue(Double.isFinite(result.aic()));
        assertTrue(result.fit().parameter("mu").effectiveDegreesOfFreedom() > 3.0);
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
