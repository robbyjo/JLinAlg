/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.junit.jupiter.api.Test;

/** Frozen base-R stats::glm and explicit delta-method reference. */
class ProbitPredictionRReferenceTest {
    @Test
    void probitFitMatchesBaseR() {
        GlmResult fit = fit();

        assertTrue(fit.converged(), fit.convergenceMessage());
        assertArrayEquals(new double[] {
            -0.051432669120702282,
            0.848776342899199121,
            -0.321522833208195247
        }, fit.coefficients(), 2e-7);
        assertArrayEquals(new double[] {
            0.54862958905233172,
            0.34821272558834554,
            0.78489110891783764
        }, fit.standardErrors(), 1e-6);
        assertEquals(13.881527713623678, fit.deviance(), 4e-10);
        assertEquals(-6.940763856811838, fit.logLikelihood(), 2e-10);
        assertEquals(19.881527713623676, fit.aic(), 4e-10);
    }

    @Test
    void expectedResponsesAndScenarioEffectsMatchBaseRDeltaMethod() {
        GlmResult fit = fit();
        ExpectedResponse[] predictions = ModelPredictions.expectedResponses(
            fit, GlmFamilies.probit(), new double[][] {
                {1, -1, 0}, {1, 0, 1}, {1, 1, 0}
            });

        assertExpected(predictions[0], -0.90020901201990144,
            0.66314328871301498, 0.18400451556218214,
            0.17641945431222442, 0.013905364160480081,
            0.65524788344827034);
        assertExpected(predictions[1], -0.37295550232889751,
            0.56793591626986184, 0.35459078118461834,
            0.21135140026419877, 0.068627732058253899,
            0.77040413576758848);
        assertExpected(predictions[2], 0.79734367377849680,
            0.63618710700993975, 0.78737426884234152,
            0.18468944026747580, 0.32651381614578667,
            0.97953544814363258);

        ScenarioDifference difference = ModelPredictions.scenarioDifference(
            fit, GlmFamilies.probit(), scenario(1.0), scenario(0.0));
        assertEquals(0.46219703457462025,
            difference.firstScenario().estimate(), 2e-7);
        assertEquals(0.53780296542537964,
            difference.secondScenario().estimate(), 2e-7);
        assertEquals(-5.5727547548854893e-6,
            difference.covarianceBetweenScenarios(), 1e-8);
        assertEquals(-0.0756059308507594, difference.difference(), 2e-7);
        assertEquals(0.18198269375813647, difference.standardError(), 1e-6);
        assertEquals(-0.43228545642628891,
            difference.confidenceLower(), 2e-6);
        assertEquals(0.28107359472477011,
            difference.confidenceUpper(), 2e-6);

        RiskRatio ratio = ModelPredictions.riskRatio(
            fit, GlmFamilies.probit(), scenario(1.0), scenario(0.0));
        assertEquals(0.85941704358034121, ratio.ratio(), 3e-7);
        assertEquals(0.31549364298201837, ratio.standardError(), 1e-6);
        assertEquals(0.36710191558183236, ratio.logStandardError(), 1e-6);
        assertEquals(0.41852966431939753, ratio.confidenceLower(), 2e-6);
        assertEquals(1.7647438587118145, ratio.confidenceUpper(), 3e-6);
    }

    @Test
    void analyticAverageMarginalEffectsMatchBaseR() {
        GlmResult fit = fit();

        AverageMarginalEffect population =
            ModelPredictions.averageMarginalEffect(
                fit, GlmFamilies.probit(), design(), 1);
        AverageMarginalEffect atAverage =
            ModelPredictions.averageMarginalEffect(
                fit, GlmFamilies.probit(), design(), null, 1,
                Averaging.AT_AVERAGE_COVARIATES, 0.95);

        assertEffect(population, 0.20008938357447134,
            0.017241281363134558, 0.16629709305540597,
            0.23388167409353672);
        assertEffect(atAverage, 0.33861276978699489,
            0.13891677881101289, 0.06634088646909281,
            0.61088465310489704);
    }

    private static void assertExpected(
            ExpectedResponse actual, double eta, double linkSe,
            double mean, double meanSe, double lower, double upper) {
        assertEquals(eta, actual.linearPredictor(), 2e-7);
        assertEquals(linkSe, actual.linkStandardError(), 1e-6);
        assertEquals(mean, actual.mean().estimate(), 2e-7);
        assertEquals(meanSe, actual.mean().standardError(), 5e-7);
        assertEquals(lower, actual.mean().confidenceLower(), 1e-6);
        assertEquals(upper, actual.mean().confidenceUpper(), 1e-6);
    }

    private static void assertEffect(
            AverageMarginalEffect actual, double estimate, double se,
            double lower, double upper) {
        assertEquals(estimate, actual.estimate(), 2e-7);
        assertEquals(se, actual.standardError(), 5e-7);
        assertEquals(lower, actual.confidenceLower(), 1e-6);
        assertEquals(upper, actual.confidenceUpper(), 1e-6);
    }

    private static GlmResult fit() {
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
