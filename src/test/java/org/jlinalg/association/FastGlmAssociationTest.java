/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
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
            assertTrue(Double.isFinite(result.pValues()[variable]));
        }
        assertEquals(2, result.parallelism());
        assertTrue(prepared.nullModel().converged());
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
