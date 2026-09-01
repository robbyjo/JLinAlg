/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class CoxRegressionTest {
    private static final double[] TIME = {1, 2, 2, 3, 4, 4, 5, 6, 7, 8};
    private static final boolean[] EVENT = {
        true, true, false, true, true, true, false, true, false, true
    };
    private static final double[][] X = {
        {-1, 0}, {-0.5, 1}, {0, 0}, {0.5, 1}, {1, 0},
        {-1, 1}, {-0.5, 0}, {0, 1}, {0.5, 0}, {1, 1}
    };

    @Test
    void matchesIndependentStatsmodelsBreslowAndEfronReferences() {
        CoxSurvivalData response = CoxSurvivalData.rightCensored(TIME, EVENT);
        CoxResult breslow = CoxRegression.fit(response, X, null,
            CoxOptions.defaults().withTies(CoxTies.BRESLOW),
            BackendPolicy.CPU);
        CoxResult efron = CoxRegression.fit(response, X, null,
            CoxOptions.defaults().withTies(CoxTies.EFRON),
            BackendPolicy.CPU);

        assertTrue(breslow.converged(), breslow.convergenceMessage());
        assertTrue(efron.converged(), efron.convergenceMessage());
        assertArrayEquals(new double[] {-0.75204494, 0.30971324},
            breslow.beta(), 1e-7);
        assertArrayEquals(new double[] {0.64103730, 0.89383259},
            breslow.standardErrors(), 1e-7);
        assertEquals(-10.209302708561657,
            breslow.logPartialLikelihood(), 1e-10);
        assertArrayEquals(new double[] {-0.79747382, 0.34327597},
            efron.beta(), 1e-7);
        assertArrayEquals(new double[] {0.64802536, 0.88900303},
            efron.standardErrors(), 1e-7);
        assertEquals(-9.952069141148542,
            efron.logPartialLikelihood(), 1e-10);
        assertEquals(Math.exp(efron.beta()[0]), efron.hazardRatios()[0], 0);
        assertTrue(efron.hazardRatioConfidenceLower()[0]
            < efron.hazardRatios()[0]);
        assertTrue(efron.hazardRatioConfidenceUpper()[0]
            > efron.hazardRatios()[0]);
    }

    @Test
    void countingProcessAndStratifiedBaselineAreConsistent() {
        int[] strata = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
        CoxSurvivalData right = CoxSurvivalData.rightCensored(
            TIME, EVENT, strata);
        CoxSurvivalData counting = new CoxSurvivalData(
            new double[TIME.length], TIME, EVENT, strata);
        CoxResult first = CoxRegression.fit(right, X, null,
            CoxOptions.defaults(), BackendPolicy.CPU);
        CoxResult second = CoxRegression.fit(counting, X, null,
            CoxOptions.defaults(), BackendPolicy.CPU);

        assertArrayEquals(first.beta(), second.beta(), 0);
        assertTrue(first.baselineHazard().stream()
            .anyMatch(point -> point.stratum() == 0));
        assertTrue(first.baselineHazard().stream()
            .anyMatch(point -> point.stratum() == 1));
        for (int index = 1; index < first.baselineHazard().size(); index++) {
            BaselineHazardPoint previous = first.baselineHazard().get(index - 1);
            BaselineHazardPoint current = first.baselineHazard().get(index);
            if (previous.stratum() == current.stratum())
                assertTrue(current.cumulativeHazard()
                    >= previous.cumulativeHazard());
        }
    }

    @Test
    void delayedEntryMatchesIndependentStatsmodelsReference() {
        // Avoid entry exactly at a failure time: statsmodels includes that
        // boundary while the counting-process convention here is (start, stop].
        double[] entry = {0, 0, 0, 0.5, 1.1, 1.2, 1.5, 2.2, 2.3, 3.2};
        CoxResult result = CoxRegression.fit(
            new CoxSurvivalData(entry, TIME, EVENT, null), X, null,
            CoxOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertArrayEquals(new double[] {-0.623897977048431, 0.484210234630982},
            result.beta(), 1e-7);
        assertArrayEquals(new double[] {0.665148745527891, 0.885114493120438},
            result.standardErrors(), 1e-7);
        assertEquals(-8.767479391027788,
            result.logPartialLikelihood(), 1e-10);
    }
}
