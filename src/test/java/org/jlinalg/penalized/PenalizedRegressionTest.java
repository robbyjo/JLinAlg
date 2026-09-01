/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.junit.jupiter.api.Test;

class PenalizedRegressionTest {
    private static final double[] RESPONSE = {2, 6, 0, 4};
    private static final double[][] ORTHOGONAL = {
        {-1, -1},
        {1, -1},
        {-1, 1},
        {1, 1}
    };

    @Test
    void ridgeMatchesOrthogonalClosedForm() {
        PenalizedRegressionResult result = PenalizedRegression.ridge(
            RESPONSE, ORTHOGONAL, 0.5);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(3.0, result.intercept(), 1e-14);
        assertArrayEquals(new double[] {4.0 / 3.0, -2.0 / 3.0},
            result.coefficients(), 2e-12);
        assertEquals(2, result.activeCoefficientCount());
        assertArrayEquals(result.fittedValues(), result.predict(ORTHOGONAL));
    }

    @Test
    void lassoMatchesOrthogonalSoftThreshold() {
        PenalizedRegressionResult result = PenalizedRegression.lasso(
            RESPONSE, ORTHOGONAL, 0.5);

        assertTrue(result.converged(), result.convergenceMessage());
        assertArrayEquals(new double[] {1.5, -0.5},
            result.coefficients(), 2e-12);
    }

    @Test
    void elasticNetMatchesOrthogonalCoordinateSolution() {
        ElasticNetOptions options = ElasticNetOptions.builder()
            .alpha(0.5)
            .standardize(false)
            .build();
        PenalizedRegressionResult result = PenalizedRegression.fit(
            RESPONSE, ORTHOGONAL, 0.5, options);

        assertTrue(result.converged(), result.convergenceMessage());
        assertArrayEquals(new double[] {1.4, -0.6},
            result.coefficients(), 2e-12);
    }

    @Test
    void ridgeConvenienceOverloadPreservesPenaltyFactors() {
        ElasticNetOptions options = ElasticNetOptions.builder()
            .alpha(0.75)
            .standardize(false)
            .penaltyFactors(1.0, 2.0)
            .build();
        PenalizedRegressionResult result = PenalizedRegression.ridge(
            RESPONSE, ORTHOGONAL, 0.5, options);

        assertEquals(0.0, result.alpha(), 0.0);
        assertArrayEquals(new double[] {4.0 / 3.0, -0.5},
            result.coefficients(), 2e-12);
    }

    @Test
    void pathUsesDescendingLambdasAndAutomaticMaximumStartsAtZero() {
        ElasticNetOptions options = ElasticNetOptions.builder()
            .alpha(1.0)
            .build();
        PenalizedRegressionPath path = PenalizedRegression.automaticPath(
            RESPONSE, ORTHOGONAL, 20, 1e-3, options);

        assertEquals(20, path.size());
        assertEquals(0, path.fit(0).activeCoefficientCount());
        assertTrue(path.fit(path.size() - 1).activeCoefficientCount() > 0);
        assertThrows(IllegalArgumentException.class,
            () -> PenalizedRegression.path(RESPONSE, ORTHOGONAL,
                new double[] {0.1, 0.2}, options));
    }

    @Test
    void crossValidationIsDeterministicAndReturnsOneSeFit() {
        int rows = 30;
        double[] response = new double[rows];
        double[][] predictors = new double[rows][2];
        for (int row = 0; row < rows; row++) {
            predictors[row][0] = row - 15.0;
            predictors[row][1] = (row % 5) - 2.0;
            response[row] = 1.5 + 2.0 * predictors[row][0]
                - 0.4 * predictors[row][1] + (row % 2 == 0 ? 0.1 : -0.1);
        }
        double[] lambdas = {2.0, 1.0, 0.5, 0.1, 0.01};
        ElasticNetOptions options = ElasticNetOptions.builder()
            .alpha(0.5)
            .build();
        PenalizedCrossValidationResult first =
            PenalizedRegressionCrossValidation.fit(
                response, predictors, lambdas, 5, 12345L, options);
        PenalizedCrossValidationResult second =
            PenalizedRegressionCrossValidation.fit(
                response, predictors, lambdas, 5, 12345L, options);

        assertArrayEquals(first.meanSquaredErrors(),
            second.meanSquaredErrors());
        assertTrue(first.lambdaOneStandardError() >= first.lambdaMinimum());
        assertEquals(first.lambdaMinimum(),
            first.minimumErrorFit().lambda(), 0.0);
        assertEquals(first.lambdaOneStandardError(),
            first.oneStandardErrorFit().lambda(), 0.0);
    }

    @Test
    void ridgeInferenceUsesEffectiveResidualDegreesOfFreedom() {
        ElasticNetOptions options = ElasticNetOptions.builder()
            .standardize(false)
            .build();
        RidgeRegressionResult result = PenalizedRegressionInference.ridge(
            RESPONSE, ORTHOGONAL, 0.5, options, BackendPolicy.CPU);

        assertArrayEquals(new double[] {4.0 / 3.0, -2.0 / 3.0},
            result.beta(), 2e-12);
        assertEquals(1.0 + 2.0 / 1.5,
            result.effectiveModelDegreesOfFreedom(), 2e-12);
        assertEquals(RESPONSE.length - result.effectiveModelDegreesOfFreedom(),
            result.residualDegreesOfFreedom(), 2e-12);
        assertEquals(DegreesOfFreedomMethod.EFFECTIVE_RESIDUAL,
            result.associationStatistics().degreesOfFreedomMethod());
        assertTrue(result.standardErrors()[0] > 0.0);
        assertTrue(result.pValues()[0] >= 0.0 && result.pValues()[0] <= 1.0);
    }

    @Test
    void activeSetRefitReturnsConditionalOlsInference() {
        ElasticNetOptions options = ElasticNetOptions.builder()
            .alpha(1.0)
            .standardize(false)
            .build();
        PenalizedRegressionResult lasso = PenalizedRegression.fit(
            RESPONSE, ORTHOGONAL, 0.5, options);
        PostSelectionOlsResult refit =
            PenalizedRegressionInference.refitActiveSet(
                RESPONSE, ORTHOGONAL, lasso, options,
                1e-12, BackendPolicy.CPU);

        assertArrayEquals(new int[] {0, 1}, refit.activePredictorIndices());
        assertArrayEquals(new double[] {3.0, 2.0, -1.0},
            refit.ols().coefficients(), 2e-12);
        assertEquals(3, refit.associationStatistics().beta().length);
    }
}
