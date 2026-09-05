/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.association.AssociationBatchResult;
import org.jlinalg.association.AssociationEngineOptions;
import org.junit.jupiter.api.Test;

class FastCoxAssociationTest {
    @Test
    void preparedScoreMatchesDirectEfficientInformation() {
        double[] time = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        boolean[] event = {
            true, false, true, true, false, true,
            true, false, true, true, false, true
        };
        double[][] covariates = new double[time.length][1];
        double[][] predictors = new double[time.length][2];
        for (int row = 0; row < time.length; row++) {
            covariates[row][0] = (row % 3) - 1.0;
            predictors[row][0] = ((row * 5) % 11) - 5.0;
            predictors[row][1] = row % 4;
        }
        CoxSurvivalData survival = CoxSurvivalData.rightCensored(time, event);
        AssociationEngineOptions options = AssociationEngineOptions.defaults()
            .withParallelism(2).withChunkSize(2);
        FastCoxAssociation prepared = FastCoxAssociation.prepare(survival,
            covariates, null, CoxOptions.defaults(), options,
            List.of("a", "a", "b", "b", "c", "c",
                "d", "d", "e", "e", "f", "f"),
            identity(time.length));

        AssociationBatchResult model = prepared.scan(predictors,
            List.of("g1", "g2"), options, CoxScoreVariance.MODEL_BASED);
        AssociationBatchResult robust = prepared.scan(predictors,
            List.of("g1", "g2"), options, CoxScoreVariance.CLUSTER_ROBUST);
        AssociationBatchResult related = prepared.scan(predictors,
            List.of("g1", "g2"), options, CoxScoreVariance.RELATEDNESS);

        assertTrue(model.successful());
        assertEquals(2, model.size());
        for (int index = 0; index < 2; index++) {
            assertTrue(Double.isFinite(model.beta()[index]));
            assertTrue(model.standardErrors()[index] > 0.0);
            assertTrue(robust.standardErrors()[index] > 0.0);
            assertTrue(related.standardErrors()[index] > 0.0);
            assertEquals(related.beta()[index], model.beta()[index], 0.0);
        }
        assertEquals(prepared.nullModel().backend(), prepared.backend());
    }

    private static double[] identity(int dimension) {
        double[] result = new double[dimension * dimension];
        for (int index = 0; index < dimension; index++)
            result[index * dimension + index] = 1.0;
        return result;
    }
}
