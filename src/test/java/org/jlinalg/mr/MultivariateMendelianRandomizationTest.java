/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class MultivariateMendelianRandomizationTest {
    @Test void recoversJointMultipleExposureMultipleOutcomeEffects() {
        MultivariateMrResult result = MultivariateMendelianRandomization.fit(
            instruments(false), List.of("x1", "x2"), List.of("y1", "y2"),
            new double[][] {{1, .55}, {.55, 1}}, false, BackendPolicy.CPU);

        assertEquals(.7, result.effect(0, 0), 2e-12);
        assertEquals(-.4, result.effect(1, 0), 2e-12);
        assertEquals(-.2, result.effect(0, 1), 2e-12);
        assertEquals(.9, result.effect(1, 1), 2e-12);
        assertEquals(4, result.overallTest().degreesOfFreedom());
        assertEquals(2, result.exposureJointTests().size());
        assertEquals(2, result.outcomeJointTests().size());
        assertEquals(0, result.cochranQ(), 1e-20);
        assertTrue(result.overallTest().pValue() < 1e-12);
        assertTrue(Math.abs(result.covariance()[2]) > 0,
            "outcome correlation must propagate into coefficient covariance");
    }

    @Test void identityOutcomeCorrelationMatchesSeparateMvmrFits() {
        List<MultivariateInstrument> joint = instruments(false);
        MultivariateMrResult result = MultivariateMendelianRandomization.fit(
            joint, List.of("x1", "x2"), List.of("y1", "y2"),
            new double[][] {{1, 0}, {0, 1}}, false, BackendPolicy.CPU);
        for (int outcome = 0; outcome < 2; outcome++) {
            List<MultivariableInstrument> separate = new ArrayList<>();
            for (MultivariateInstrument value : joint)
                separate.add(new MultivariableInstrument(value.variantId(),
                    value.exposureEffects(), value.exposureStandardErrors(),
                    value.outcomeEffects()[outcome],
                    value.outcomeStandardErrors()[outcome]));
            MultivariableMrResult expected =
                MultivariableMendelianRandomization.fit(separate,
                    List.of("x1", "x2"), false, BackendPolicy.CPU);
            assertArrayEquals(expected.beta(), new double[] {
                result.effect(0, outcome), result.effect(1, outcome)}, 2e-12);
        }
    }

    @Test void multiplicativeRandomEffectsInflateUncertainHeterogeneousFits() {
        List<MultivariateInstrument> values = instruments(true);
        MultivariateMrResult fixed = MultivariateMendelianRandomization.fit(
            values, List.of("x1", "x2"), List.of("y1", "y2"),
            new double[][] {{1, .25}, {.25, 1}}, false, BackendPolicy.CPU);
        MultivariateMrResult random = MultivariateMendelianRandomization.fit(
            values, List.of("x1", "x2"), List.of("y1", "y2"),
            new double[][] {{1, .25}, {.25, 1}}, true, BackendPolicy.CPU);
        assertTrue(fixed.heterogeneityPValue() < .05);
        assertTrue(random.dispersion() > 1);
        assertTrue(random.standardErrors()[0] > fixed.standardErrors()[0]);
    }

    @Test void rejectsMvmrAndInvalidCorrelationShapesAtTheBoundary() {
        List<MultivariateInstrument> values = instruments(false);
        assertThrows(IllegalArgumentException.class,
            () -> MultivariateMendelianRandomization.ivw(values,
                List.of("x1", "x2"), List.of("only-one"),
                new double[][] {{1}}));
        assertThrows(IllegalArgumentException.class,
            () -> MultivariateMendelianRandomization.ivw(values,
                List.of("x1", "x2"), List.of("y1", "y2"),
                new double[][] {{1, 1.1}, {1.1, 1}}));
        assertThrows(IllegalArgumentException.class,
            () -> MultivariateMendelianRandomization.ivw(values,
                List.of("x1", "x2"), List.of("y1", "y2"),
                new double[][] {{1, 0}, null}));
    }

    @Test void multivariatePressoFindsAJointOutcomeOutlierReproducibly() {
        List<MultivariateInstrument> values = instruments(true);
        MultivariateMrPressoResult first = MultivariateMrPresso.analyze(values,
            List.of("x1", "x2"), List.of("y1", "y2"),
            new double[][] {{1, .25}, {.25, 1}}, .2, 59, 73L,
            BackendPolicy.CPU);
        MultivariateMrPressoResult second = MultivariateMrPresso.analyze(values,
            List.of("x1", "x2"), List.of("y1", "y2"),
            new double[][] {{1, .25}, {.25, 1}}, .2, 59, 73L,
            BackendPolicy.CPU);
        assertEquals(first.globalPValue(), second.globalPValue());
        assertTrue(first.globalPValue() > 0,
            "finite Monte Carlo tests must not report an impossible zero p-value");
        assertTrue(first.instrumentTests().stream()
            .allMatch(test -> test.empiricalPValue() > 0));
        assertTrue(first.outlierVariants().contains("v1"));
        assertNotNull(first.outlierCorrectedEstimate());
        assertEquals(.7, first.outlierCorrectedEstimate().effect(0, 0), 2e-12);
        assertEquals(.9, first.outlierCorrectedEstimate().effect(1, 1), 2e-12);
    }

    private static List<MultivariateInstrument> instruments(boolean outlier) {
        double[][] x = {
            {.10, .02}, {.02, .12}, {.15, -.03}, {-.04, .16},
            {.12, .10}, {-.08, .07}, {.18, -.06}, {.04, -.11}
        };
        List<MultivariateInstrument> result = new ArrayList<>();
        for (int i = 0; i < x.length; i++) {
            double y1 = .7 * x[i][0] - .4 * x[i][1];
            double y2 = -.2 * x[i][0] + .9 * x[i][1];
            if (outlier && i == 0) { y1 += .3; y2 -= .25; }
            result.add(new MultivariateInstrument("v" + (i + 1), x[i],
                new double[] {.01, .012}, new double[] {y1, y2},
                new double[] {.02 + .001 * i, .025 + .0005 * i}));
        }
        return result;
    }
}
