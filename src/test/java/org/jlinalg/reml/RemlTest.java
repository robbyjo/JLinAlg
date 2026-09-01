/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.junit.jupiter.api.Test;

class RemlTest {
    @Test
    void identityOnlyRemlMatchesOlsResidualVariance() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[][] fixed = {
            {1.0, 0.0},
            {1.0, 1.0},
            {1.0, 2.0},
            {1.0, 3.0},
            {1.0, 4.0}
        };

        RemlResult result = Reml.fit(response, fixed,
            List.of(VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder().initialVariances(1.0).build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.7 / 3.0, result.varianceComponents()[0], 1e-7);
        assertArrayEquals(new double[] {0.6, 2.1},
            result.fixedEffects(), 1e-10);
        assertEquals("cpu", result.backend().selectedBackend());
        assertEquals(2.0, result.degreesOfFreedom()[0], 0.0);
        assertEquals(2, result.beta().length);
        assertEquals(2, result.tStatistics().length);
        assertEquals(2, result.pValues().length);
    }

    @Test
    void profileMlUsesMaximumLikelihoodResidualVariance() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[][] fixed = {
            {1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0},
            {1.0, 3.0}, {1.0, 4.0}
        };

        RemlResult result = Reml.fit(response, fixed,
            List.of(VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder().initialVariances(1.0)
                .varianceEstimation(VarianceEstimation.ML).build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.7 / 5.0, result.varianceComponents()[0], 1e-7);
        assertEquals(VarianceEstimation.ML, result.varianceEstimation());
        assertEquals(result.logLikelihood(), result.restrictedLogLikelihood());
    }

    @Test
    void satterthwaiteMatchesExactResidualDfForIdentityCovariance() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[][] fixed = {
            {1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0},
            {1.0, 3.0}, {1.0, 4.0}
        };

        RemlResult result = Reml.fit(response, fixed,
            List.of(VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder()
                .initialVariances(1.0)
                .degreesOfFreedomMethod(DegreesOfFreedomMethod.SATTERTHWAITE)
                .build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertArrayEquals(new double[] {3.0, 3.0},
            result.degreesOfFreedom(), 2e-8);
        assertEquals(DegreesOfFreedomMethod.SATTERTHWAITE,
            result.associationStatistics().degreesOfFreedomMethod());
    }

    @Test
    void kenwardRogerMatchesExactLinearModelInferenceForIdentityCovariance() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[][] fixed = {
            {1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0},
            {1.0, 3.0}, {1.0, 4.0}
        };

        RemlResult result = Reml.fit(response, fixed,
            List.of(VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder()
                .initialVariances(1.0)
                .degreesOfFreedomMethod(DegreesOfFreedomMethod.KENWARD_ROGER)
                .build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertArrayEquals(new double[] {3.0, 3.0},
            result.degreesOfFreedom(), 2e-8);
        assertArrayEquals(result.fixedEffectCovariance(),
            result.fixedEffectInferenceCovariance(), 2e-10);
        assertEquals(DegreesOfFreedomMethod.KENWARD_ROGER,
            result.associationStatistics().degreesOfFreedomMethod());
    }

    @Test
    void kenwardRogerInflatesCovarianceForUnbalancedRandomIntercept() {
        double[] response = {1.2, 2.1, 4.2, 3.9, 5.1, 7.8, 8.4, 9.1};
        double[][] fixed = {
            {1, 0}, {1, 1},
            {1, 0}, {1, 1}, {1, 2},
            {1, 0}, {1, 1}, {1, 2}
        };
        List<String> groups = List.of(
            "a", "a", "b", "b", "b", "c", "c", "c");

        RemlResult result = Reml.fit(response, fixed,
            List.of(
                VarianceComponent.randomIntercept("group", groups),
                VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder()
                .initialVariances(4.0, 1.0)
                .degreesOfFreedomMethod(DegreesOfFreedomMethod.KENWARD_ROGER)
                .build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        double[] model = result.fixedEffectCovariance();
        double[] adjusted = result.fixedEffectInferenceCovariance();
        assertTrue(adjusted[0] > model[0]);
        assertTrue(adjusted[3] > model[3]);
        assertEquals(Math.sqrt(adjusted[0]), result.standardErrors()[0], 1e-12);
        assertTrue(Double.isFinite(result.degreesOfFreedom()[0]));
        assertTrue(result.degreesOfFreedom()[0] > 0.0);
    }

    @Test
    void randomEffectFactoriesBuildExpectedCovariances() {
        List<String> groups = List.of("a", "a", "b");
        assertArrayEquals(new double[] {
            1, 1, 0,
            1, 1, 0,
            0, 0, 1
        }, VarianceComponent.randomIntercept("group", groups).covariance());
        assertArrayEquals(new double[] {
            1, 2, 0,
            2, 4, 0,
            0, 0, 9
        }, VarianceComponent.randomSlope(
            "slope", groups, new double[] {1, 2, 3}).covariance());
        assertArrayEquals(new double[] {
            1, 1, 0,
            1, 1, 0,
            0, 0, 1
        }, VarianceComponent.fromRandomEffectDesign("design", new double[][] {
            {1, 0}, {1, 0}, {0, 1}
        }).covariance());
    }

    @Test
    void balancedRandomInterceptMatchesClosedFormReml() {
        double[] response = {0.0, 1.0, 2.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0};
        double[][] fixed = new double[response.length][1];
        for (double[] row : fixed) {
            row[0] = 1.0;
        }
        double[] groupCovariance = new double[response.length * response.length];
        for (int row = 0; row < response.length; row++) {
            for (int column = 0; column < response.length; column++) {
                if (row / 3 == column / 3) {
                    groupCovariance[row * response.length + column] = 1.0;
                }
            }
        }

        RemlResult result = Reml.fit(response, fixed,
            List.of(
                new VarianceComponent("group", response.length, groupCovariance),
                VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder().initialVariances(10.0, 2.0).build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(47.0 / 3.0, result.varianceComponents()[0], 1e-5);
        assertEquals(1.0, result.varianceComponents()[1], 1e-6);
        assertEquals(5.0, result.fixedEffects()[0], 1e-10);
        assertEquals(2, result.componentNames().size());
    }

    @Test
    void rejectsNonsymmetricCovarianceComponent() {
        assertThrows(IllegalArgumentException.class,
            () -> new VarianceComponent("bad", 2,
                new double[] {1.0, 1.0, 0.0, 1.0}));
    }

    @Test
    void supportsKnownCovarianceContribution() {
        double[] response = {1.0, 2.0, 5.0, 7.0, 9.0};
        double[] design = {
            1.0, 0.0, 1.0, 1.0, 1.0, 2.0,
            1.0, 3.0, 1.0, 4.0
        };
        double[] known = new double[response.length * response.length];
        for (int index = 0; index < response.length; index++) {
            known[index * response.length + index] = 0.1;
        }

        RemlResult result = Reml.fitWithKnownCovariance(
            response, design, response.length, 2,
            List.of(VarianceComponent.identity("estimated", response.length)),
            known,
            RemlOptions.builder().initialVariances(0.2).build(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.7 / 3.0 - 0.1,
            result.varianceComponents()[0], 1e-7);
    }
}
