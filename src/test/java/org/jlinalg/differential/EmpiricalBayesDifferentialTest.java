/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

final class EmpiricalBayesDifferentialTest {
    private static final double[][] DESIGN = {
        {1, 0}, {1, 0}, {1, 0}, {1, 1}, {1, 1}, {1, 1}
    };
    private static final double[] CONTRAST = {0, 1};

    @Test void continuousModerationSharesInformationAcrossFeatures() {
        double[][] data = {
            {1.0, 1.2, 0.8, 3.0, 3.1, 2.9},
            {2.0, 2.1, 1.9, 2.1, 2.0, 2.2},
            {4.0, 4.3, 3.7, 5.0, 5.2, 4.8},
            {8.0, 7.0, 9.0, 8.1, 7.2, 8.9}
        };
        DifferentialFit fit = EmpiricalBayesDifferential.fitContinuous(
            data, DESIGN, CONTRAST);
        assertEquals("limma", fit.method());
        assertEquals(2.0, fit.results().get(0).effect(), 1e-12);
        assertTrue(fit.results().get(0).pValue() < 1e-4);
        assertTrue(fit.priorDegreesOfFreedom() > 0.0);
        assertTrue(fit.results().stream().allMatch(result ->
            result.moderatedVariance() > 0.0
                && Double.isFinite(result.standardError())));
    }

    @Test void voomAndNegativeBinomialKeepSeparateEstimatorContracts() {
        double[][] counts = {
            {10, 12, 11, 30, 33, 29},
            {100, 110, 90, 105, 95, 100},
            {4, 5, 3, 16, 14, 15},
            {40, 35, 45, 42, 38, 41}
        };
        DifferentialFit voom = EmpiricalBayesDifferential.fitVoom(
            counts, DESIGN, CONTRAST);
        double[][] normalizedFixture = {
            {10, 12, 11, 30, 33, 29},
            {100, 100, 100, 100, 100, 100},
            {200, 200, 200, 200, 200, 200},
            {50, 50, 50, 50, 50, 50},
            {80, 80, 80, 80, 80, 80}
        };
        DifferentialFit negativeBinomial =
            EmpiricalBayesDifferential.fitNegativeBinomial(
                normalizedFixture, DESIGN, CONTRAST);
        assertEquals("voom", voom.method());
        assertEquals("negative-binomial", negativeBinomial.method());
        assertTrue(voom.results().stream().allMatch(result ->
            result.meanPrecisionWeight() > 0.0));
        assertEquals(Math.log(30.666666666666668 / 11.0),
            negativeBinomial.results().get(0).effect(), 0.12);
        assertTrue(negativeBinomial.results().stream().allMatch(result ->
            result.dispersion() > 0.0 && result.converged()));
        double logGeometricMean = 0.0;
        for (double factor : negativeBinomial.sizeFactors())
            logGeometricMean += Math.log(factor) / negativeBinomial.sizeFactors().length;
        assertEquals(0.0, logGeometricMean, 1e-12);
    }
}
