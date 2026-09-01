/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.ols;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;
import org.junit.jupiter.api.Test;

class OlsTest {
    @Test
    void fitsKnownFullRankRegression() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[][] design = {
            {1.0, 0.0},
            {1.0, 1.0},
            {1.0, 2.0},
            {1.0, 3.0},
            {1.0, 4.0}
        };

        OlsResult result = Ols.fit(response, design,
            OlsOptions.defaults(), BackendPolicy.CPU);

        assertArrayEquals(new double[] {0.6, 2.1},
            result.coefficients(), 1e-12);
        assertEquals(0.7, result.residualSumOfSquares(), 1e-12);
        assertEquals(0.7 / 3.0, result.residualVariance(), 1e-12);
        assertEquals(2, result.rank());
        assertEquals(3, result.residualDegreesOfFreedom());
        assertFalse(result.rankDeficient());
        assertFalse(result.minimumNorm());

        double[] covariance = result.covariance();
        assertEquals(0.14, covariance[0], 1e-12);
        assertEquals(-0.7 / 15.0, covariance[1], 1e-12);
        assertEquals(-0.7 / 15.0, covariance[2], 1e-12);
        assertEquals(0.7 / 30.0, covariance[3], 1e-12);
        assertEquals("cpu", result.backend().selectedBackend());
        assertTrue(result.pValues()[1] < 0.001);
        assertEquals(result.tStatistics()[1] * result.tStatistics()[1],
            result.testContrast(new double[][] {{0.0, 1.0}}).statistic(), 1e-10);
    }

    @Test
    void rejectsRankDeficiencyByDefault() {
        double[] response = {2.0, 4.0, 6.0};
        double[][] design = {
            {1.0, 1.0},
            {1.0, 1.0},
            {1.0, 1.0}
        };

        IllegalArgumentException exception = assertThrows(
            IllegalArgumentException.class,
            () -> Ols.fit(response, design,
                OlsOptions.defaults(), BackendPolicy.CPU));
        assertTrue(exception.getMessage().contains("rank deficient"));
    }

    @Test
    void computesMinimumNormSolutionWhenRequested() {
        double[] response = {2.0, 4.0, 6.0};
        double[][] design = {
            {1.0, 1.0},
            {1.0, 1.0},
            {1.0, 1.0}
        };
        OlsOptions options = new OlsOptions(
            RankDeficiencyStrategy.MINIMUM_NORM, 0.95);

        OlsResult result = Ols.fit(
            response, design, options, BackendPolicy.CPU);

        assertArrayEquals(new double[] {2.0, 2.0},
            result.coefficients(), 1e-12);
        assertEquals(1, result.rank());
        assertEquals(2, result.residualDegreesOfFreedom());
        assertTrue(result.rankDeficient());
        assertTrue(result.minimumNorm());
        assertArrayEquals(new double[] {4.0, 4.0, 4.0},
            result.fittedValues(), 1e-12);
    }

    @Test
    void resultArraysAreDefensiveCopies() {
        OlsResult result = Ols.fit(
            new double[] {1.0, 2.0, 4.0},
            new double[][] {{1.0}, {1.0}, {1.0}},
            OlsOptions.defaults(), BackendPolicy.CPU);
        double original = result.coefficients()[0];
        double[] exposed = result.coefficients();
        exposed[0] = -100.0;
        assertEquals(original, result.coefficients()[0]);
    }

    @Test
    void supportsWeightsOffsetAndCompleteCaseOmission() {
        double[] response = {11.0, 13.0, Double.NaN, 17.0};
        double[][] design = {
            {1.0, 0.0}, {1.0, 1.0},
            {1.0, Double.NaN}, {1.0, 3.0}
        };
        double[] weights = {1.0, 2.0, 3.0, 4.0};
        double[] offset = {10.0, 10.0, 10.0, 10.0};
        OlsOptions options = new OlsOptions(
            RankDeficiencyStrategy.ERROR, 0.95, MissingDataPolicy.OMIT);

        OlsResult result = Ols.fit(response, design, weights, offset,
            options, BackendPolicy.CPU);

        assertArrayEquals(new double[] {1.0, 2.0},
            result.coefficients(), 1e-12);
        assertEquals(3, result.observations());
        assertEquals(4, result.originalObservations());
        assertEquals(1, result.omittedObservations());
        assertArrayEquals(new int[] {0, 1, 3}, result.retainedRows());
    }
}
