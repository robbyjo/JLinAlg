/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class ProbitGlmTest {
    @Test
    void retainsNormalTailInformationWhenReportedProbabilityRounds() {
        GlmFamily family = GlmFamilies.binomialProbit();

        assertSame(family, GlmFamilies.probit());
        assertEquals(7.619853024160527e-24,
            family.inverseLink(-10.0), 2e-38);
        assertEquals(1.0, family.inverseLink(10.0));
        assertEquals(7.69459862670642e-23,
            family.meanDerivative(10.0), 2e-37);
        assertEquals(7.619853024160533e-24,
            family.varianceAtPredictor(10.0, 1.0), 2e-37);
        assertEquals(7.770077433040107e-22,
            family.workingWeight(0.0, 10.0, 1.0, 1.0), 2e-35);
        assertEquals(-53.23128515051247,
            family.logLikelihoodAtPredictor(0.0, 10.0, 1.0, 1.0, 1.0),
            2e-13);
        assertEquals(106.46257030102494,
            family.unitDevianceAtPredictor(0.0, 10.0, 1.0), 4e-13);

        assertEquals(5.725571222524577e-300,
            family.inverseLink(-37.0), 4e-313);
        assertTrue(family.varianceAtPredictor(37.0, 1.0) > 0.0);
        assertTrue(family.workingWeight(0.0, 37.0, 1.0, 1.0) > 0.0);
    }

    @Test
    void fitsGroupedBinomialProportionsAndRejectsInvalidResponses() {
        double[] response = {0.1, 0.3, 0.6, 0.8};
        double[] trials = {10, 10, 10, 10};
        double[][] design = {{1, -1}, {1, 0}, {1, 1}, {1, 2}};

        GlmResult fit = Glm.fit(response, design, GlmFamilies.probit(),
            trials, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(fit.converged(), fit.convergenceMessage());
        assertTrue(Double.isFinite(fit.logLikelihood()));
        assertEquals("binomial(probit)", fit.family());
        assertThrows(IllegalArgumentException.class,
            () -> Glm.fit(new double[] {0.0, 1.01},
                new double[][] {{1.0}, {1.0}}, GlmFamilies.probit()));
    }
}
