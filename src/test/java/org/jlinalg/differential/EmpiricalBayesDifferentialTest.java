/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

final class EmpiricalBayesDifferentialTest {
    @Test void exactCountFitsUseUnscaledInformationAndFailedFitsSuppressInference() {
        var fit = EmpiricalBayesDifferential.fitNegativeBinomial(new double[][] {
            {10,10,10,10,10,10}, {20,20,20,20,20,20},
            {20,20,20,40,40,40}, {0,0,0,0,0,0}}, DESIGN, CONTRAST);
        var changed = fit.results().get(2);
        // Base-R inverse X'WX at the fitted dispersion, not residual-scaled OLS.
        assertEquals(.23990279696049335, changed.standardError(), 1e-12);
        assertEquals(.0038612085215580121, changed.pValue(), 1e-13);
        var zero = fit.results().get(3);
        assertTrue(!zero.converged());
        assertTrue(Double.isNaN(zero.standardError()));
        assertTrue(Double.isNaN(zero.statistic()));
        assertTrue(Double.isNaN(zero.pValue()));
    }

    @Test void voomAcceptsTwoFeatures() {
        var fit = EmpiricalBayesDifferential.fitVoom(new double[][] {
            {10,20,15,25}, {20,30,25,35}},
            new double[][] {{1,0},{1,0},{1,1},{1,1}}, CONTRAST);
        assertEquals(2, fit.results().size());
        assertTrue(fit.results().stream().allMatch(r -> Double.isFinite(r.pValue())));
    }

    @Test void exactGaussianFeatureStillReceivesModeratedPriorUncertainty() {
        var fit = EmpiricalBayesDifferential.fitContinuous(new double[][] {
            {1,1,1,3,3,3}, {2,3,1,3,4,2}, {4,6,2,7,9,5}}, DESIGN, CONTRAST);
        var row = fit.results().get(0);
        assertEquals(Math.sqrt((2.0/3) * row.moderatedVariance()), row.standardError(), 1e-12);
        assertTrue(row.standardError() > 1e-6);
    }
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
