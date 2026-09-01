/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceEstimation;
import org.junit.jupiter.api.Test;

class SparseLinearMixedModelTest {
    @Test
    void sparseRandomInterceptMatchesDenseReferenceFit() {
        double[] response = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        double[][] fixed = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        RandomEffectTerm subject = RandomEffectTerm.randomIntercept(
            "subject", List.of("a", "a", "a", "b", "b", "b",
                "c", "c", "c"));
        RemlOptions options = RemlOptions.builder()
            .initialVariances(10.0, 2.0)
            .maximumIterations(200).build();

        LinearMixedModelResult dense = LinearMixedModel.fit(response, fixed,
            List.of(subject), options, BackendPolicy.CPU);
        SparseLinearMixedModelResult sparse = SparseLinearMixedModel.fit(
            response, fixed, List.of(subject), options, BackendPolicy.CPU);

        assertArrayEquals(dense.reml().varianceComponents(),
            sparse.varianceComponents(), 2e-5);
        assertArrayEquals(dense.beta(), sparse.beta(), 1e-8);
        assertArrayEquals(dense.standardErrors(), sparse.standardErrors(), 2e-6);
        assertArrayEquals(dense.randomEffects("subject").estimates(),
            sparse.randomEffects("subject").estimates(), 2e-5);
        assertEquals(dense.reml().restrictedLogLikelihood(),
            sparse.logLikelihood(), 1e-7);
        assertTrue(sparse.converged());
        assertEquals(3, sparse.randomCoefficientCount());
        assertTrue(sparse.equationNonzeroCount() <= 3);
        assertTrue(sparse.factorNonzeroCount() <= 3);
    }

    @Test
    void supportsCrossedIndependentSparseTermsWithoutDenseCovariance() {
        double[] response = {2, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 14};
        double[][] fixed = {
            {1, 0}, {1, 1}, {1, 2}, {1, 3},
            {1, 0}, {1, 1}, {1, 2}, {1, 3},
            {1, 0}, {1, 1}, {1, 2}, {1, 3}
        };
        RandomEffectTerm subject = RandomEffectTerm.randomIntercept("subject",
            List.of("a", "a", "a", "a", "b", "b", "b", "b",
                "c", "c", "c", "c"));
        RandomEffectTerm batch = RandomEffectTerm.randomIntercept("batch",
            List.of("x", "y", "x", "y", "x", "y", "x", "y",
                "x", "y", "x", "y"));
        RemlOptions options = RemlOptions.builder()
            .initialVariances(2.0, 1.0, 1.0)
            .maximumIterations(200).build();

        SparseLinearMixedModelResult result = SparseLinearMixedModel.fit(
            response, fixed, List.of(subject, batch), options,
            BackendPolicy.CPU);

        assertEquals(3, result.varianceComponents().length);
        assertEquals(5, result.randomCoefficientCount());
        assertEquals(2, result.beta().length);
        assertTrue(Double.isFinite(result.logLikelihood()));
        assertTrue(result.equationNonzeroCount() < 25);
    }

    @Test
    void comparesNestedMaximumLikelihoodModels() {
        double[] response = {1.0, 2.2, 3.1, 4.4, 2.0, 3.1, 4.3, 5.5,
            0.5, 1.8, 2.9, 4.0};
        double[][] reducedFixed = {
            {1}, {1}, {1}, {1}, {1}, {1},
            {1}, {1}, {1}, {1}, {1}, {1}
        };
        double[][] fullFixed = {
            {1, 0}, {1, 1}, {1, 2}, {1, 3},
            {1, 0}, {1, 1}, {1, 2}, {1, 3},
            {1, 0}, {1, 1}, {1, 2}, {1, 3}
        };
        RandomEffectTerm subject = RandomEffectTerm.randomIntercept("subject",
            List.of("a", "a", "a", "a", "b", "b", "b", "b",
                "c", "c", "c", "c"));
        RemlOptions options = RemlOptions.builder()
            .initialVariances(1, 1)
            .varianceEstimation(VarianceEstimation.ML).build();
        SparseLinearMixedModelResult reduced = SparseLinearMixedModel.fit(
            response, reducedFixed, List.of(subject), options,
            BackendPolicy.CPU);
        SparseLinearMixedModelResult full = SparseLinearMixedModel.fit(
            response, fullFixed, List.of(subject), options,
            BackendPolicy.CPU);

        MixedModelComparisonResult comparison =
            MixedModelComparison.compare(reduced, full);
        assertEquals(1, comparison.degreesOfFreedom());
        assertTrue(comparison.likelihoodRatioStatistic() >= 0.0);
        assertTrue(comparison.pValue() >= 0.0
            && comparison.pValue() <= 1.0);
    }
}
