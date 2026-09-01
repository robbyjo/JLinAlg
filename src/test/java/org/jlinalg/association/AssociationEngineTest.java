/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.PValueScale;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.junit.jupiter.api.Test;

class AssociationEngineTest {
    private static final AssociationEngineOptions PARALLEL =
        new AssociationEngineOptions(3, 1, BackendPolicy.CPU,
            AssociationFailurePolicy.FAIL_FAST,
            VariableMissingPolicy.MEAN_IMPUTE);

    @Test
    void genericPredictorScanMatchesIndependentOlsFitsAndKeepsOrder() {
        double[] response = {2.0, 3.5, 5.2, 6.7, 7.9, 10.1, 11.2, 13.4};
        double[][] covariates = {
            {1, 0}, {1, 1}, {1, 0}, {1, 1},
            {1, 0}, {1, 1}, {1, 0}, {1, 1}
        };
        double[][] predictors = {
            {0, 2, 0}, {1, 1, 0}, {2, 0, 1}, {3, 2, 1},
            {4, 1, 2}, {5, 0, 3}, {6, 2, 5}, {7, 1, 8}
        };

        AssociationBatchResult result = ParallelAssociationEngine.scanPredictors(
            response, covariates, predictors, List.of("linear", "cycle", "fib"),
            AssociationModels.ols(OlsOptions.defaults()), PARALLEL);

        assertEquals(List.of("linear", "cycle", "fib"), result.names());
        assertEquals(3, result.parallelism());
        assertArrayEquals(result.beta(), result.effectSizes(), 0.0);
        assertEquals(Math.log10(result.pValues()[0]),
            result.pValues(PValueScale.LOG10)[0], 1e-13);
        assertEquals(-Math.log10(result.pValues()[0]),
            result.pValues(PValueScale.NEGATIVE_LOG10)[0], 1e-13);
        assertEquals(result.beta()[0], result.estimate(0).beta(), 0.0);
        for (int predictor = 0; predictor < 3; predictor++) {
            double[][] full = append(covariates, predictors, predictor);
            var expected = Ols.fit(response, full,
                OlsOptions.defaults(), BackendPolicy.CPU);
            assertEquals(expected.coefficients()[2], result.beta()[predictor], 1e-11);
            assertEquals(expected.standardErrors()[2],
                result.standardErrors()[predictor], 1e-11);
            assertEquals(expected.tStatistics()[2],
                result.statistics()[predictor], 1e-11);
            assertEquals(expected.pValues()[2], result.pValues()[predictor], 1e-11);
        }
    }

    @Test
    void genericResponseScanMatchesIndependentFits() {
        double[][] design = {
            {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}
        };
        double[][] responses = {
            {1.1, 7.0, 2.0}, {2.8, 6.1, 2.1}, {5.2, 5.1, 3.0},
            {6.9, 3.8, 4.2}, {9.1, 3.0, 4.8}, {10.9, 2.2, 6.1}
        };

        AssociationBatchResult result = ParallelAssociationEngine.scanResponses(
            responses, design, 1, List.of("a", "b", "c"),
            AssociationModels.ols(OlsOptions.defaults()), PARALLEL);

        for (int trait = 0; trait < 3; trait++) {
            double[] response = column(responses, trait);
            var expected = Ols.fit(response, design,
                OlsOptions.defaults(), BackendPolicy.CPU);
            assertEquals(expected.coefficients()[1], result.beta()[trait], 1e-12);
            assertEquals(expected.standardErrors()[1],
                result.standardErrors()[trait], 1e-12);
        }
    }

    @Test
    void recordNanContinuesAfterNonEstimablePredictor() {
        double[] response = {1, 2, 4, 7, 11};
        double[][] covariates = {{1}, {1}, {1}, {1}, {1}};
        double[][] predictors = {
            {1, 0}, {1, 1}, {1, 2}, {1, 4}, {1, 8}
        };
        AssociationEngineOptions options = new AssociationEngineOptions(
            2, 1, BackendPolicy.CPU, AssociationFailurePolicy.RECORD_NAN,
            VariableMissingPolicy.ERROR);

        AssociationBatchResult result = ParallelAssociationEngine.scanPredictors(
            response, covariates, predictors, List.of("constant", "usable"),
            AssociationModels.ols(OlsOptions.defaults()), options);

        assertTrue(Double.isNaN(result.beta()[0]));
        assertTrue(Double.isFinite(result.beta()[1]));
        assertEquals(1, result.failures().size());
        assertEquals("constant", result.failures().get(0).name());

        double[][] withAllMissing = {
            {Double.NaN, 0}, {Double.NaN, 1}, {Double.NaN, 2},
            {Double.NaN, 4}, {Double.NaN, 8}
        };
        AssociationBatchResult fast = FastOlsAssociation.scanPredictors(
            response, covariates, withAllMissing, List.of("missing", "usable"),
            null, null, OlsOptions.defaults(),
            options.withChunkSize(2));
        assertTrue(Double.isNaN(fast.beta()[0]));
        assertTrue(Double.isFinite(fast.beta()[1]));
        assertEquals("missing", fast.failures().get(0).name());
    }

    @Test
    void fastOlsPredictorAndResponseScansMatchFullWeightedFits() {
        double[] response = {12.1, 12.9, 15.4, 16.8, 19.1, 20.7, 23.2, 24.8};
        double[][] covariates = {
            {1, 0}, {1, 1}, {1, 0}, {1, 1},
            {1, 0}, {1, 1}, {1, 0}, {1, 1}
        };
        double[][] predictors = {
            {0, 0}, {1, 1}, {2, Double.NaN}, {3, 1},
            {4, 2}, {5, 3}, {6, 5}, {7, 8}
        };
        double[] weights = {1, 2, 1, 3, 2, 1, 2, 4};
        double[] offset = {10, 10, 10, 10, 10, 10, 10, 10};

        AssociationBatchResult fast = FastOlsAssociation.scanPredictors(
            response, covariates, predictors, List.of("x", "z"), weights,
            offset, OlsOptions.defaults(), PARALLEL);

        for (int predictor = 0; predictor < 2; predictor++) {
            double[][] imputed = imputedAppend(covariates, predictors, predictor);
            var expected = Ols.fit(response, imputed, weights, offset,
                OlsOptions.defaults(), BackendPolicy.CPU);
            assertEquals(expected.coefficients()[2], fast.beta()[predictor], 1e-10);
            assertEquals(expected.standardErrors()[2],
                fast.standardErrors()[predictor], 1e-10);
            assertEquals(expected.pValues()[2], fast.pValues()[predictor], 1e-10);
        }

        double[][] traits = {
            {12.1, 3.1}, {12.9, 4.0}, {15.4, 5.2}, {16.8, 6.1},
            {19.1, 8.2}, {20.7, 9.1}, {23.2, 11.0}, {24.8, 11.8}
        };
        AssociationBatchResult responseScan = FastOlsAssociation.scanResponses(
            traits, covariates, 1, List.of("first", "second"), weights,
            null, OlsOptions.defaults(), PARALLEL);
        for (int trait = 0; trait < 2; trait++) {
            var expected = Ols.fit(column(traits, trait), covariates,
                weights, null, OlsOptions.defaults(), BackendPolicy.CPU);
            assertEquals(expected.coefficients()[1],
                responseScan.beta()[trait], 1e-11);
            assertEquals(expected.standardErrors()[1],
                responseScan.standardErrors()[trait], 1e-11);
        }
        assertArrayEquals(new double[] {6, 6},
            responseScan.degreesOfFreedom(), 0.0);
    }

    private static double[][] append(
            double[][] fixed, double[][] predictors, int predictor) {
        double[][] result = new double[fixed.length][fixed[0].length + 1];
        for (int row = 0; row < fixed.length; row++) {
            System.arraycopy(fixed[row], 0, result[row], 0, fixed[row].length);
            result[row][fixed[row].length] = predictors[row][predictor];
        }
        return result;
    }

    private static double[][] imputedAppend(
            double[][] fixed, double[][] predictors, int predictor) {
        double sum = 0;
        int count = 0;
        for (double[] row : predictors) {
            if (Double.isFinite(row[predictor])) {
                sum += row[predictor];
                count++;
            }
        }
        double[][] result = append(fixed, predictors, predictor);
        for (double[] row : result) {
            if (!Double.isFinite(row[fixed[0].length]))
                row[fixed[0].length] = sum / count;
        }
        return result;
    }

    private static double[] column(double[][] matrix, int column) {
        double[] result = new double[matrix.length];
        for (int row = 0; row < matrix.length; row++)
            result[row] = matrix[row][column];
        return result;
    }
}
