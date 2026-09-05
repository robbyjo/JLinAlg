/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import jdistlib.Normal;
import org.junit.jupiter.api.Test;

class StatisticalTestsTest {
    private static final double TOLERANCE = 1e-9;

    @Test
    void fisherExactMatchesCanonicalRReference() {
        StatisticalTestResult result = StatisticalTests.fisherExact(new long[][] {
            {1, 9},
            {11, 3}
        });

        assertEquals(0.0027594561852200836, result.pValue(), 1e-14);
        assertEquals("exact conditional minimum-likelihood", result.pValueMethod());
        assertTrue(result.estimates().get("odds ratio") > 0.0);
        assertTrue(result.confidenceInterval().orElseThrow().upper() < 1.0);

        StatisticalTestResult rectangular = StatisticalTests.fisherExact(new long[][] {
            {1, 2, 1},
            {2, 1, 2}
        });
        assertTrue(rectangular.pValue() > 0.0 && rectangular.pValue() <= 1.0);
        assertTrue(rectangular.parameters().get("enumerated tables") > 1.0);
    }

    @Test
    void blockedRankTestsUseBaseRStatistics() {
        double[][] blocks = {
            {1, 2, 3},
            {2, 1, 3},
            {1, 3, 2},
            {2, 3, 1}
        };

        StatisticalTestResult friedman = StatisticalTests.friedman(blocks);
        StatisticalTestResult quade = StatisticalTests.quade(blocks);

        assertEquals(1.5, friedman.statistic(), TOLERANCE);
        assertEquals(Math.exp(-0.75), friedman.pValue(), TOLERANCE);
        assertEquals(0.6923076923076923, quade.statistic(), TOLERANCE);
        assertEquals(2.0, quade.parameters().get("num df"), 0.0);
        assertEquals(6.0, quade.parameters().get("denom df"), 0.0);
    }

    @Test
    void correlationSupportsPearsonAndExactRankTests() {
        double[] increasing = {1, 2, 3, 4, 5};
        double[] decreasing = {5, 4, 3, 2, 1};
        StatisticalTestResult spearman = StatisticalTests.correlation(
            increasing, decreasing, CorrelationMethod.SPEARMAN,
            Alternative.TWO_SIDED);
        StatisticalTestResult kendall = StatisticalTests.correlation(
            increasing, decreasing, CorrelationMethod.KENDALL,
            Alternative.TWO_SIDED);
        StatisticalTestResult pearson = StatisticalTests.correlation(
            increasing, new double[] {2, 1, 4, 3, 5});

        assertEquals(-1.0, spearman.estimates().get("rho"), TOLERANCE);
        assertEquals(1.0 / 60.0, spearman.pValue(), TOLERANCE);
        assertEquals(-1.0, kendall.estimates().get("tau"), TOLERANCE);
        assertEquals(1.0 / 60.0, kendall.pValue(), TOLERANCE);
        assertEquals(0.8, pearson.estimates().get("correlation"), TOLERANCE);
        assertEquals(2.3094010767585034, pearson.statistic(), TOLERANCE);
        assertTrue(pearson.confidenceInterval().isPresent());
    }

    @Test
    void oneWayAnovaSupportsClassicalAndWelchTests() {
        double[][] groups = {
            {1, 2, 3},
            {2, 4, 6},
            {5, 6, 7}
        };

        StatisticalTestResult classical = StatisticalTests.oneWayAnova(groups, true);
        StatisticalTestResult welch = StatisticalTests.oneWayAnova(groups);

        assertEquals(6.0, classical.statistic(), TOLERANCE);
        assertEquals(1.0 / 27.0, classical.pValue(), TOLERANCE);
        assertEquals(10.204724409448819, welch.statistic(), TOLERANCE);
        assertEquals(3.789473684210526, welch.parameters().get("denom df"), TOLERANCE);
    }

    @Test
    void proportionAndTrendMatchBaseRFormulas() {
        StatisticalTestResult oneSample = StatisticalTests.proportion(
            new int[] {5}, new int[] {10}, null, Alternative.TWO_SIDED,
            0.95, true);
        StatisticalTestResult trend = StatisticalTests.proportionTrend(
            new int[] {10, 20, 30}, new int[] {100, 100, 100});

        assertEquals(0.0, oneSample.statistic(), TOLERANCE);
        assertEquals(1.0, oneSample.parameters().get("df"), 0.0);
        assertEquals(1.0, oneSample.pValue(), TOLERANCE);
        assertEquals(0.5, oneSample.estimates().get("proportion"), TOLERANCE);
        assertTrue(oneSample.confidenceInterval().orElseThrow().lower() > 0.2);
        assertTrue(oneSample.confidenceInterval().orElseThrow().upper() < 0.8);
        assertEquals(12.5, trend.statistic(), TOLERANCE);
        assertEquals(ChiSquareReference.upperTailOneDf(12.5), trend.pValue(), 1e-12);
    }

    @Test
    void mcnemarAndMantelHaenszelReturnNamedUniformResults() {
        StatisticalTestResult mcnemar = StatisticalTests.mcnemar(new long[][] {
            {10, 20},
            {30, 40}
        });
        StatisticalTestResult mantel = StatisticalTests.mantelHaenszel(new long[][][] {
            {{12, 5}, {7, 10}},
            {{8, 9}, {4, 13}}
        });

        assertEquals(1.62, mcnemar.statistic(), TOLERANCE);
        assertEquals(1.0, mcnemar.parameters().get("df"), 0.0);
        assertEquals("Mantel-Haenszel X-squared", mantel.statisticName());
        assertTrue(mantel.estimates().get("common odds ratio") > 1.0);
        assertTrue(mantel.confidenceInterval().isPresent());

        StatisticalTestResult general = StatisticalTests.mantelHaenszel(new long[][][] {
            {{10, 5, 3}, {4, 8, 6}},
            {{7, 6, 4}, {3, 9, 5}}
        });
        assertEquals("Cochran-Mantel-Haenszel M^2", general.statisticName());
        assertEquals(2.0, general.parameters().get("df"), 0.0);
    }

    @Test
    void anovaPowerMatchesDocumentedBaseRExampleAndSolvesN() {
        AnovaPowerResult direct = StatisticalTests.powerAnova(4, 5, 1, 3, 0.05);
        AnovaPowerResult solved = StatisticalTests.powerAnova(
            4.0, null, 1.0, 3.0, 0.05, 0.80);

        assertEquals(0.3535594, direct.power(), 5e-7);
        assertEquals(11.92613, solved.observationsPerGroup(), 5e-5);
    }

    @Test
    void jdistlibBackedFamiliesShareTheUniformContract() {
        double[] first = {1, 2, 3, 4, 6};
        double[] second = {2, 3, 5, 7, 9};
        double[] values = {1, 2, 3, 4, 6, 2, 3, 5, 7, 9};
        int[] groups = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

        StatisticalTestResult[] results = {
            StatisticalTests.ansariBradley(first, second),
            StatisticalTests.bartlett(values, groups),
            StatisticalTests.binomial(7, 10, 0.5),
            StatisticalTests.chiSquareGoodnessOfFit(
                new long[] {10, 20, 30}, new double[] {1, 1, 1}),
            StatisticalTests.chiSquareIndependence(new long[][] {{10, 20}, {20, 10}}),
            StatisticalTests.flignerKilleen(values, groups),
            StatisticalTests.kolmogorovSmirnov(first, second),
            StatisticalTests.kolmogorovSmirnov(
                first, new Normal(), Alternative.TWO_SIDED, true),
            StatisticalTests.kruskalWallis(values, groups),
            StatisticalTests.mood(first, second),
            StatisticalTests.poisson(7, 2.0, 3.0),
            StatisticalTests.poisson(7, 5, 2.0, 3.0, 1.0, Alternative.TWO_SIDED),
            StatisticalTests.shapiroWilk(first),
            StatisticalTests.studentT(first, 0.0),
            StatisticalTests.studentT(first, second, 0.0, false,
                Alternative.TWO_SIDED, 0.95),
            StatisticalTests.pairedT(first, second, 0.0,
                Alternative.TWO_SIDED, 0.95),
            StatisticalTests.variance(first, second),
            StatisticalTests.wilcoxonSignedRank(
                first, 0.0, Alternative.TWO_SIDED, true),
            StatisticalTests.wilcoxonRankSum(
                first, second, 0.0, Alternative.TWO_SIDED, true)
        };

        for (StatisticalTestResult result : results) {
            assertFalse(result.method().isBlank());
            assertFalse(result.statisticName().isBlank());
            assertTrue(Double.isNaN(result.pValue())
                || result.pValue() >= 0.0 && result.pValue() <= 1.0);
            assertFalse(result.pValueMethod().isBlank());
        }
    }

    @Test
    void invalidAndNonFiniteInputsAreRejected() {
        assertThrows(IllegalArgumentException.class,
            () -> StatisticalTests.studentT(new double[] {1, Double.NaN}, 0));
        assertThrows(IllegalArgumentException.class,
            () -> StatisticalTests.fisherExact(new long[][] {{1, -2}, {4, 5}}));
        assertThrows(IllegalArgumentException.class,
            () -> StatisticalTests.friedman(new double[][] {{1, 2}, {3}}));
    }

    /** Independent identity for a chi-square(1) upper tail. */
    private static final class ChiSquareReference {
        private ChiSquareReference() { }

        static double upperTailOneDf(double value) {
            return 2.0 * Normal.cumulative(Math.sqrt(value), 0.0, 1.0, false, false);
        }
    }
}
