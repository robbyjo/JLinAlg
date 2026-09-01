/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;
import org.junit.jupiter.api.Test;

class GlmTest {
    @Test
    void gaussianIdentityMatchesOls() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[][] design = {
            {1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0},
            {1.0, 3.0}, {1.0, 4.0}
        };

        GlmResult result = Glm.fit(response, design, GlmFamilies.gaussian(),
            null, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertArrayEquals(new double[] {0.6, 2.1},
            result.coefficients(), 1e-12);
        assertEquals(0.7, result.deviance(), 1e-12);
        assertEquals(0.7 / 3.0, result.dispersion(), 1e-12);
    }

    @Test
    void binomialInterceptMatchesObservedProportion() {
        double[] response = {0.0, 0.0, 0.0, 1.0, 1.0};
        double[][] intercept = {{1.0}, {1.0}, {1.0}, {1.0}, {1.0}};

        GlmResult result = Glm.fit(response, intercept, GlmFamilies.binomial(),
            null, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(Math.log(0.4 / 0.6), result.coefficients()[0], 1e-10);
        for (double fitted : result.fittedMeans()) {
            assertEquals(0.4, fitted, 1e-10);
        }
        assertEquals(1.0, result.dispersion());
    }

    @Test
    void poissonInterceptMatchesSampleMean() {
        double[] response = {0.0, 1.0, 2.0, 3.0};
        double[][] intercept = {{1.0}, {1.0}, {1.0}, {1.0}};

        GlmResult result = Glm.fit(response, intercept, GlmFamilies.poisson(),
            null, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(Math.log(1.5), result.coefficients()[0], 1e-10);
        for (double fitted : result.fittedMeans()) {
            assertEquals(1.5, fitted, 1e-10);
        }
    }

    @Test
    void binomialRejectsOutOfRangeResponse() {
        assertThrows(IllegalArgumentException.class,
            () -> Glm.fit(new double[] {0.0, 2.0},
                new double[][] {{1.0}, {1.0}}, GlmFamilies.binomial(),
                null, null, GlmOptions.defaults(), BackendPolicy.CPU));
    }

    @Test
    void poissonOffsetRepresentsExposure() {
        double[] response = {1.0, 2.0, 2.0, 4.0};
        double[] exposure = {1.0, 2.0, 2.0, 4.0};
        double[] offset = new double[exposure.length];
        for (int index = 0; index < exposure.length; index++) {
            offset[index] = Math.log(exposure[index]);
        }
        double[][] intercept = {{1.0}, {1.0}, {1.0}, {1.0}};

        GlmResult result = Glm.fit(response, intercept, GlmFamilies.poisson(),
            null, offset, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.0, result.coefficients()[0], 1e-10);
        assertArrayEquals(exposure, result.fittedMeans(), 1e-10);
    }

    @Test
    void gammaLogInterceptMatchesLogMeanForEqualWeights() {
        double[] response = {1.0, 2.0, 3.0, 6.0};
        double[][] intercept = {{1.0}, {1.0}, {1.0}, {1.0}};

        GlmResult result = Glm.fit(response, intercept, GlmFamilies.gamma(),
            null, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(Math.log(3.0), result.coefficients()[0], 1e-9);
    }

    @Test
    void negativeBinomialInterceptMatchesSampleMean() {
        double[] response = {0.0, 1.0, 2.0, 5.0};
        double[][] intercept = {{1.0}, {1.0}, {1.0}, {1.0}};

        GlmResult result = Glm.fit(response, intercept,
            GlmFamilies.negativeBinomial(2.0), null, null,
            GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(Math.log(2.0), result.coefficients()[0], 1e-9);
    }

    @Test
    void omitsNonFiniteRowsWhenRequested() {
        double[] response = {0.0, 1.0, Double.NaN, 1.0};
        double[][] intercept = {{1.0}, {1.0}, {Double.NaN}, {1.0}};
        GlmOptions options = GlmOptions.builder()
            .missingDataPolicy(MissingDataPolicy.OMIT).build();

        GlmResult result = Glm.fit(response, intercept, GlmFamilies.binomial(),
            null, null, options, BackendPolicy.CPU);

        assertEquals(3, result.observations());
        assertEquals(1, result.omittedObservations());
        assertArrayEquals(new int[] {0, 1, 3}, result.retainedRows());
    }
}
