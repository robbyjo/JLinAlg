/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gee.Gee;
import org.jlinalg.gee.GeeOptions;
import org.jlinalg.gee.GeeResult;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.ols.RankDeficiencyStrategy;
import org.junit.jupiter.api.Test;

class ModelPredictionsTest {
    @Test
    void populationAverageAndAverageCovariateAreDistinctEstimands() {
        GlmResult fit = probitFit();
        double[][] first = scenario(1.0);
        double[][] second = scenario(0.0);

        ScenarioDifference population = ModelPredictions.scenarioDifference(
            fit, GlmFamilies.probit(), first, second, null, null,
            Averaging.POPULATION_AVERAGE, 0.95);
        ScenarioDifference atAverage = ModelPredictions.scenarioDifference(
            fit, GlmFamilies.probit(), first, second, null, null,
            Averaging.AT_AVERAGE_COVARIATES, 0.95);

        assertNotEquals(population.difference(), atAverage.difference());
        assertEquals(-0.0756059308507594, population.difference(), 2e-7);
        assertEquals(-0.12771868501149303, atAverage.difference(), 2e-7);
        assertTrue(Double.isFinite(population.covarianceBetweenScenarios()));
    }

    @Test
    void riskRatioRejectsNonpositiveResponseMeans() {
        double[] response = {-2, -1, 0, 1, 2};
        double[][] design = {{1, -2}, {1, -1}, {1, 0}, {1, 1}, {1, 2}};
        GlmResult fit = Glm.fit(response, design, GlmFamilies.gaussian(),
            null, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.riskRatio(fit, GlmFamilies.gaussian(),
                new double[][] {{1, -2}}, new double[][] {{1, 2}}));
    }

    @Test
    void geeUsesTheSameExpectedResponseContract() {
        double[] response = {0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0};
        double[][] design = new double[response.length][2];
        int[] cluster = new int[response.length];
        for (int row = 0; row < response.length; row++) {
            design[row][0] = 1.0;
            design[row][1] = row % 2;
            cluster[row] = row / 2;
        }
        GeeResult fit = Gee.fit(response, design, cluster, null,
            GlmFamilies.binomial(), null, null,
            GeeOptions.defaults(), BackendPolicy.CPU);

        ExpectedResponse[] predictions = ModelPredictions.expectedResponses(
            fit, GlmFamilies.binomial(),
            new double[][] {{1, 0}, {1, 1}});

        assertEquals(2, predictions.length);
        assertTrue(predictions[0].mean().estimate() > 0.0);
        assertTrue(predictions[0].mean().estimate() < 1.0);
        assertTrue(predictions[0].mean().confidenceLower()
            <= predictions[0].mean().estimate());
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.expectedResponses(fit,
                GlmFamilies.probit(), new double[][] {{1, 0}}));
    }

    @Test
    void rejectsMismatchedFamiliesAndUnconvergedFits() {
        GlmResult fit = probitFit();
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.expectedResponses(fit,
                GlmFamilies.binomial(), new double[][] {{1, 0, 0}}));

        GlmResult unfinished = Glm.fit(
            new double[] {1, 5, 3, 8},
            new double[][] {{1}, {1}, {1}, {1}},
            GlmFamilies.poisson(), null, null,
            GlmOptions.builder().maximumIterations(1)
                .initialCoefficients(-2).build(), BackendPolicy.CPU);
        assertTrue(!unfinished.converged());
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.expectedResponses(unfinished,
                GlmFamilies.poisson(), new double[][] {{1}}));

        double[] response = {0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0};
        double[][] design = new double[response.length][2];
        int[] cluster = new int[response.length];
        for (int row = 0; row < response.length; row++) {
            design[row][0] = 1.0;
            design[row][1] = row % 2;
            cluster[row] = row / 2;
        }
        GeeResult unfinishedGee = Gee.fit(response, design, cluster, null,
            GlmFamilies.binomial(), null, null,
            GeeOptions.builder().maximumIterations(1)
                .relativeTolerance(1e-30)
                .initialCoefficients(4, -4).build(), BackendPolicy.CPU);
        assertTrue(!unfinishedGee.converged());
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.expectedResponses(unfinishedGee,
                GlmFamilies.binomial(), new double[][] {{1, 0}}));
    }

    @Test
    void rankDeficientGlmRejectsNonEstimablePredictionTargets() {
        double[][] design = {
            {1, 0, 0}, {1, 1, 1}, {1, 2, 2},
            {1, 3, 3}, {1, 4, 4}
        };
        GlmResult fit = Glm.fit(new double[] {1, 2, 5, 7, 9}, design,
            GlmFamilies.gaussian(), null, null,
            GlmOptions.builder().rankDeficiencyStrategy(
                RankDeficiencyStrategy.MINIMUM_NORM).build(),
            BackendPolicy.CPU);

        fit.requireEstimable(new double[] {1, 1, 1});
        assertEquals(1, ModelPredictions.expectedResponses(fit,
            GlmFamilies.gaussian(), new double[][] {{1, 1, 1}}).length);
        assertThrows(IllegalArgumentException.class,
            () -> fit.requireEstimable(new double[] {1, 1, 0}));
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.expectedResponses(fit,
                GlmFamilies.gaussian(), new double[][] {{1, 1, 0}}));
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.scenarioDifference(fit,
                GlmFamilies.gaussian(), new double[][] {{1, 1, 1}},
                new double[][] {{1, 1, 0}}));
        assertThrows(IllegalArgumentException.class,
            () -> ModelPredictions.averageMarginalEffect(fit,
                GlmFamilies.gaussian(), design, 1));
    }

    private static GlmResult probitFit() {
        return Glm.fit(RESPONSE, design(), GlmFamilies.probit(),
            null, null, GlmOptions.builder()
                .relativeTolerance(1e-12).build(), BackendPolicy.CPU);
    }

    private static double[][] design() {
        double[][] result = new double[X.length][3];
        for (int row = 0; row < result.length; row++) {
            result[row] = new double[] {1.0, X[row], row % 2};
        }
        return result;
    }

    private static double[][] scenario(double z) {
        double[][] result = design();
        for (double[] row : result) row[2] = z;
        return result;
    }

    private static final double[] X = {
        -2, -1.7, -1.4, -1.1, -.8, -.5, -.2, .1,
        .4, .7, 1, 1.3, 1.6, 1.9, 2.2, 2.5
    };
    private static final double[] RESPONSE = {
        0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1
    };
}
