/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import jdistlib.Normal;
import jdistlib.T;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.settest.SetTestScoreState;
import org.junit.jupiter.api.Test;

class FastGlmAssociationTest {
    @Test
    void binomialScoreScanMatchesDirectEfficientScoreCalculation() {
        double[] response = {0, 0, 0, 1, 1, 1, 0, 1};
        double[][] covariates = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        double[][] predictors = {
            {0, 2}, {0, 1}, {1, 0}, {1, 2},
            {2, 1}, {2, 0}, {0, 1}, {2, 2}
        };
        AssociationEngineOptions execution =
            new AssociationEngineOptions(2, 1, BackendPolicy.CPU,
                AssociationFailurePolicy.FAIL_FAST,
                VariableMissingPolicy.MEAN_IMPUTE);

        FastGlmAssociation prepared = FastGlmAssociation.prepare(
            response, covariates, GlmFamilies.binomial(), null, null,
            GlmOptions.defaults(), execution);
        AssociationBatchResult result = prepared.scan(
            predictors, List.of("g1", "g2"), execution);

        for (int variable = 0; variable < 2; variable++) {
            double mean = 0.0;
            for (double[] row : predictors) mean += row[variable];
            mean /= predictors.length;
            double score = 0.0;
            double information = 0.0;
            for (int row = 0; row < response.length; row++) {
                double centered = predictors[row][variable] - mean;
                score += centered * (response[row] - 0.5);
                information += 0.25 * centered * centered;
            }
            double expectedBeta = score / information;
            double expectedSe = Math.sqrt(1.0 / information);
            assertEquals(expectedBeta, result.effectSizes()[variable], 1e-12);
            assertEquals(expectedSe, result.standardErrors()[variable], 1e-12);
            assertEquals(expectedBeta / expectedSe,
                result.tOrZStatistics()[variable], 1e-12);
            double expectedP = 2.0 * Normal.cumulative(
                Math.abs(expectedBeta / expectedSe),
                0.0, 1.0, false, false);
            assertEquals(expectedP, result.pValues()[variable], 1e-15);
            assertEquals(Double.POSITIVE_INFINITY,
                result.degreesOfFreedom()[variable]);
        }
        assertEquals(2, result.parallelism());
        assertTrue(prepared.nullModel().converged());
    }

    @Test
    void gaussianScoreScanUsesResidualStudentTInference() {
        double[] response = {1, 2, 2, 4, 5, 7, 7, 9};
        double[][] covariates = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        double[][] predictor = {
            {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}
        };
        AssociationEngineOptions execution =
            AssociationEngineOptions.defaults()
                .withBackendPolicy(BackendPolicy.CPU);

        AssociationBatchResult result = FastGlmAssociation.prepare(
            response, covariates, GlmFamilies.gaussian(), null, null,
            GlmOptions.defaults(), execution).scan(
                predictor, List.of("x"), execution);

        double statistic = result.tOrZStatistics()[0];
        double degreesOfFreedom = response.length - 2.0;
        double expectedP = 2.0 * T.cumulative(
            Math.abs(statistic), degreesOfFreedom, false, false);
        assertEquals(degreesOfFreedom, result.degreesOfFreedom()[0], 0.0);
        assertEquals(expectedP, result.pValues()[0], 1e-15);
    }

    @Test
    void probitScorePreparationRetainsRoundedExtremeTailInformation() {
        double[] response = {1, 1, 1, 1, 0, 0, 0, 0};
        double[][] covariates = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        double[] offset = {10, 10, 10, 10, -10, -10, -10, -10};
        double[][] predictor = {
            {0}, {1}, {2}, {3}, {0}, {1}, {2}, {3}
        };
        AssociationEngineOptions execution =
            AssociationEngineOptions.defaults()
                .withBackendPolicy(BackendPolicy.CPU);

        FastGlmAssociation prepared = FastGlmAssociation.prepare(
            response, covariates, GlmFamilies.probit(), null, offset,
            GlmOptions.defaults(), execution);
        assertEquals(1.0, prepared.nullModel().fittedMeans()[0], 0.0);

        AssociationBatchResult result = prepared.scan(
            predictor, List.of("x"), execution);
        assertTrue(result.failures().isEmpty());
        assertTrue(Double.isFinite(result.effectSizes()[0]));
        assertTrue(Double.isFinite(result.standardErrors()[0]));
        assertTrue(Double.isFinite(result.pValues()[0]));

        SetTestScoreState score = prepared.score(new double[][] {
            {0, 1, 2, 3, 0, 1, 2, 3}
        });
        assertTrue(Double.isFinite(score.scores()[0]));
        assertTrue(Double.isFinite(score.information()[0]));
        assertTrue(score.information()[0] > 0.0);
    }

    @Test
    void recordsInvalidPredictorAndContinuesWithinSameBlock() {
        double[] response = {0, 0, 1, 1, 0, 1};
        double[][] covariates = {{1}, {1}, {1}, {1}, {1}, {1}};
        double[][] predictors = {
            {Double.NaN, 0}, {Double.NaN, 1}, {Double.NaN, 0},
            {Double.NaN, 1}, {Double.NaN, 2}, {Double.NaN, 2}
        };
        AssociationEngineOptions execution =
            new AssociationEngineOptions(1, 2, BackendPolicy.CPU,
                AssociationFailurePolicy.RECORD_NAN,
                VariableMissingPolicy.MEAN_IMPUTE);
        FastGlmAssociation prepared = FastGlmAssociation.prepare(
            response, covariates, GlmFamilies.binomial(), null, null,
            GlmOptions.defaults(), execution);

        AssociationBatchResult result = prepared.scan(
            predictors, List.of("missing", "usable"), execution);

        assertTrue(Double.isNaN(result.beta()[0]));
        assertTrue(Double.isFinite(result.beta()[1]));
        assertEquals(1, result.failures().size());
        assertEquals("missing", result.failures().get(0).name());
    }
}
