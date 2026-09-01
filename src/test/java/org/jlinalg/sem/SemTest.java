/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class SemTest {
    @Test
    void covarianceMlRecoversObservedRegressionPath() {
        SemModel model = SemModel.builder("x", "y")
            .regression("y", "x", 0.3)
            .variance("x", 1.0)
            .variance("y", 1.0)
            .build();
        double[] covariance = {
            1.0, 0.6,
            0.6, 1.36
        };
        SemFitResult result = Sem.fitCovariance(covariance, 1_000, model,
            SemOptions.defaults(), BackendPolicy.CPU);
        assertTrue(result.converged());
        assertEquals(0.6, result.parameter("y~x").estimate(), 1e-5);
        assertEquals(1.0, result.parameter("x~~x").estimate(), 1e-5);
        assertEquals(1.0, result.parameter("y~~y").estimate(), 1e-5);
        assertTrue(result.parameter("y~x").standardError() > 0.0);
        assertEquals(0.0, result.srmr(), 1e-5);
    }

    @Test
    void pathModelReportsPositiveDegreesOfFreedomAndFitIndices() {
        // x variance 1; m=.5x+e_m(.75); y=.7m+e_y(.51)
        double[] covariance = {
            1.0, 0.5, 0.35,
            0.5, 1.0, 0.7,
            0.35, 0.7, 1.0
        };
        SemModel model = SemModel.builder("x", "m", "y")
            .regression("m", "x", 0.4)
            .regression("y", "m", 0.6)
            .variance("x", 1.0)
            .variance("m", 0.8)
            .variance("y", 0.5)
            .build();
        SemFitResult result = Sem.fitCovariance(covariance, 2_000, model,
            SemOptions.defaults(), BackendPolicy.CPU);
        assertEquals(1, result.degreesOfFreedom());
        assertEquals(0.5, result.parameter("m~x").estimate(), 1e-4);
        assertEquals(0.7, result.parameter("y~m").estimate(), 1e-4);
        assertTrue(result.cfi() > 0.99);
        assertTrue(result.rmsea() < 0.01);
    }

    @Test
    void repeatedLabelsCreateEqualityConstraint() {
        SemModel model = SemModel.builder("x", "m", "y")
            .regression("equal", "m", "x", 0.4)
            .regression("equal", "y", "m", 0.4)
            .build();
        assertEquals(4, model.freeParameterCount());
    }
}
