/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class ArimaErrorLinearMixedModelTest {
    @Test
    void profileRemlRecoversArOneResidualCorrelation() {
        int observations = 36;
        List<String> groups = new ArrayList<>(observations);
        double[] groupEffects = {-1.0, 0.2, 0.8};
        double[] response = new double[observations];
        double[][] fixed = new double[observations][1];
        Random random = new Random(424242L);
        double error = 0.0;
        for (int time = 0; time < observations; time++) {
            int group = time % groupEffects.length;
            groups.add("g" + group);
            fixed[time][0] = 1.0;
            error = 0.55 * error + random.nextGaussian() * 0.6;
            response[time] = 2.0 + groupEffects[group] + error;
        }
        ArimaErrorLmmOptions options = ArimaErrorLmmOptions.builder()
            .maximumFunctionEvaluations(100)
            .optimizationTolerance(2e-3)
            .remlOptions(RemlOptions.builder()
                .initialVariances(1.0, 0.5)
                .build())
            .build();

        ArimaErrorLmmResult result = ArimaErrorLinearMixedModel.fit(
            response, fixed,
            List.of(RandomEffectTerm.randomIntercept("group", groups)),
            ArimaOrder.ar(1), options, BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.55, result.autoregressive()[0], 0.25);
        assertEquals(2.0, result.beta()[0], 0.5);
        assertEquals(0, result.differencingLoss());
        assertTrue(result.standardErrors()[0] > 0.0);
    }

    @Test
    void integratedErrorsDifferenceEveryModelMatrix() {
        int observations = 30;
        List<String> groups = new ArrayList<>(observations);
        double[] response = new double[observations];
        double[][] trend = new double[observations][1];
        double[] effects = {-0.6, 0.6};
        Random random = new Random(1234L);
        double walk = 0.0;
        for (int time = 0; time < observations; time++) {
            int group = time % 2;
            groups.add("g" + group);
            trend[time][0] = time;
            walk += random.nextGaussian() * 0.2;
            response[time] = 3.0 + 1.75 * time + effects[group] + walk;
        }
        ArimaErrorLmmOptions options = ArimaErrorLmmOptions.builder()
            .remlOptions(RemlOptions.builder()
                .initialVariances(0.5, 0.2)
                .build())
            .build();

        ArimaErrorLmmResult result = ArimaErrorLinearMixedModel.fit(
            response, trend,
            List.of(RandomEffectTerm.randomIntercept("group", groups)),
            ArimaOrder.arima(0, 1, 0), options, BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertTrue(result.differenced());
        assertEquals(1, result.differencingLoss());
        assertEquals(1.75, result.beta()[0], 0.15);
    }
}
