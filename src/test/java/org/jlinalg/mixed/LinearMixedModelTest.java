/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class LinearMixedModelTest {
    @Test
    void balancedRandomInterceptReturnsConditionalModesAndInference() {
        double[] response = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        double[][] fixed = new double[response.length][1];
        for (double[] row : fixed) {
            row[0] = 1.0;
        }
        List<String> groups = List.of(
            "a", "a", "a", "b", "b", "b", "c", "c", "c");
        RandomEffectTerm subject = RandomEffectTerm.randomIntercept(
            "subject", groups);

        LinearMixedModelResult result = LinearMixedModel.fit(
            response, fixed, List.of(subject),
            RemlOptions.builder().initialVariances(10.0, 2.0).build(),
            BackendPolicy.CPU);

        assertTrue(result.reml().converged(),
            result.reml().convergenceMessage());
        assertArrayEquals(new double[] {47.0 / 3.0, 1.0},
            result.reml().varianceComponents(), 1e-5);
        assertArrayEquals(new double[] {-47.0 / 12.0, 0.0, 47.0 / 12.0},
            result.randomEffects("subject").estimates(), 2e-6);
        assertEquals(3, result.randomEffects("subject")
            .predictionErrorVariances().length);
        assertEquals(response.length - 1 - 1,
            result.reml().degreesOfFreedom()[0]);
        assertEquals(1, result.associationStatistics().pValues().length);
        assertTrue(Math.abs(result.conditionalResiduals()[0])
            < Math.abs(result.reml().residuals()[0]));
        assertArrayEquals(result.beta(), result.fixef(), 0.0);
        assertEquals(result.randomEffects("subject"),
            result.ranef().get("subject"));
        assertEquals(2, result.varCorr().size());
        assertArrayEquals(result.conditionalFittedValues(),
            result.fittedValues(), 0.0);
    }

    @Test
    void groupedSlopeUsesOneCoefficientPerGroup() {
        RandomEffectTerm slope = RandomEffectTerm.randomSlope(
            "site-age", List.of("a", "a", "b"),
            new double[] {1.0, 2.0, 3.0});

        assertEquals(2, slope.coefficients());
        assertEquals(List.of("a", "b"), slope.coefficientNames());
        assertTrue(slope.sparse());
        assertEquals(3, slope.nonzeroCount());
        assertArrayEquals(new int[] {0, 1, 2, 3}, slope.rowPointers());
        assertArrayEquals(new int[] {0, 0, 1}, slope.columnIndices());
        assertArrayEquals(new double[] {1, 2, 3}, slope.sparseValues());
        assertArrayEquals(new double[] {
            1, 0,
            2, 0,
            0, 3
        }, slope.design());
    }

    @Test
    void preparedRefitAndNewLevelPredictionMatchLme4Semantics() {
        double[] response = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        double[][] fixed = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        RandomEffectTerm subject = RandomEffectTerm.randomIntercept("subject",
            List.of("a", "a", "a", "b", "b", "b", "c", "c", "c"));
        PreparedLinearMixedModel prepared = new PreparedLinearMixedModel(
            fixed, List.of(subject), RemlOptions.builder()
                .initialVariances(10, 2).build(), BackendPolicy.CPU);
        LinearMixedModelResult original = prepared.fit(response);
        double[] shifted = java.util.Arrays.stream(response)
            .map(value -> value + 2.0).toArray();
        LinearMixedModelResult refitted = prepared.refit(original, shifted);
        assertEquals(original.beta()[0] + 2.0, refitted.beta()[0], 1e-7);

        double[][] newFixed = {{1}, {1}};
        RandomEffectTerm known = RandomEffectTerm.randomIntercept(
            "subject", List.of("a", "new"));
        double[] marginal = MixedModelPrediction.marginal(original, newFixed);
        double[] conditional = MixedModelPrediction.conditional(
            original, newFixed, List.of(known), true);
        assertEquals(marginal[1], conditional[1], 1e-12);
        assertEquals(marginal[0] + original.randomEffects("subject").estimates()[0],
            conditional[0], 1e-12);
        org.junit.jupiter.api.Assertions.assertThrows(IllegalArgumentException.class,
            () -> MixedModelPrediction.conditional(
                original, newFixed, List.of(known), false));
    }
}
