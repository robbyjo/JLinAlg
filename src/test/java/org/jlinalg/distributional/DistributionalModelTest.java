/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PSplineTerm;
import org.jlinalg.gam.PenalizedPredictor;
import org.junit.jupiter.api.Test;

class DistributionalModelTest {
    @Test
    void gaussianLocationScaleModelsMeanAndHeteroscedasticity() {
        int observations = 240;
        double[] response = new double[observations];
        double[][] design = new double[observations][2];
        for (int row = 0; row < observations; row++) {
            double x = -1.0 + 2.0 * row / (observations - 1.0);
            double sigma = Math.exp(-0.15 + 0.45 * x);
            double error = Math.sqrt(2.0) * Math.sin(19.0 * row + 0.3);
            design[row][0] = 1.0;
            design[row][1] = x;
            response[row] = 1.25 + 1.8 * x + sigma * error;
        }

        DistributionalResult result = DistributionalModel.fit(response,
            List.of(PenalizedPredictor.linear(design),
                PenalizedPredictor.linear(design)),
            DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(1.25, result.parameter("mu").coefficients()[0], 0.03);
        assertEquals(1.8, result.parameter("mu").coefficients()[1], 0.16);
        assertEquals(-0.15, result.parameter("sigma").coefficients()[0], 0.04);
        assertEquals(0.45, result.parameter("sigma").coefficients()[1], 0.06);
    }

    @Test
    void multinomialUsesOnePredictorPerNonbaselineCategory() {
        double[] response = new double[100];
        double[][] intercept = new double[100][1];
        for (int row = 0; row < response.length; row++) {
            intercept[row][0] = 1.0;
            response[row] = row < 50 ? 0.0 : row < 80 ? 1.0 : 2.0;
        }
        PenalizedPredictor predictor = PenalizedPredictor.linear(intercept);

        DistributionalResult result = DistributionalModel.fit(response,
            List.of(predictor, predictor),
            DistributionalFamilies.multinomial(3),
            DistributionalOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(Math.log(50.0 / 20.0),
            result.parameter("logit0").coefficients()[0], 1e-8);
        assertEquals(Math.log(30.0 / 20.0),
            result.parameter("logit1").coefficients()[0], 1e-8);
    }

    @Test
    void locationAndScaleCanEachUsePenalizedSmooths() {
        int observations = 100;
        double[] x = new double[observations];
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            x[row] = row / (observations - 1.0);
            intercept[row][0] = 1.0;
            double mean = Math.sin(2.0 * Math.PI * x[row]);
            double sigma = Math.exp(-1.0
                + 0.4 * Math.cos(2.0 * Math.PI * x[row]));
            response[row] = mean + sigma * Math.sqrt(2.0)
                * Math.sin(23.0 * row + 0.2);
        }
        PenalizedPredictor mean = PenalizedPredictor.additive(intercept,
            List.of(PSplineTerm.of("s.mean(x)", x, 8)),
            new double[] {1.0}, BackendPolicy.CPU);
        PenalizedPredictor scale = PenalizedPredictor.additive(intercept,
            List.of(PSplineTerm.of("s.scale(x)", x, 7)),
            new double[] {2.0}, BackendPolicy.CPU);

        DistributionalResult result = DistributionalModel.fit(response,
            List.of(mean, scale),
            DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertTrue(result.parameter("mu").effectiveDegreesOfFreedom() > 3.0);
        assertTrue(result.parameter("sigma").effectiveDegreesOfFreedom() > 2.0);
        assertTrue(result.logLikelihood() > -80.0);
    }
}
