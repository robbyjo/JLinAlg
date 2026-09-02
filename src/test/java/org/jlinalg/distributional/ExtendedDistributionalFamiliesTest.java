/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;
import org.junit.jupiter.api.Test;

class ExtendedDistributionalFamiliesTest {
    @Test
    void specialFunctionsMatchKnownValues() {
        assertEquals(Math.log(24.0), SpecialFunctions.logGamma(5.0), 1e-12);
        assertEquals(-0.5772156649015329,
            SpecialFunctions.digamma(1.0), 2e-10);
        assertEquals(Math.PI * Math.PI / 6.0,
            SpecialFunctions.trigamma(1.0), 2e-10);
    }

    @Test
    void gammaMeanAndShapeFitJointly() {
        int observations = 180;
        double[] response = new double[observations];
        double[][] meanDesign = new double[observations][2];
        double[][] shapeIntercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            double x = -1.0 + 2.0 * row / (observations - 1.0);
            double mean = Math.exp(0.4 + 0.3 * x);
            double noise = Math.exp(0.32 * Math.sqrt(2.0)
                * Math.sin(17.0 * row + 0.4) - 0.32 * 0.32 / 2.0);
            meanDesign[row][0] = 1.0;
            meanDesign[row][1] = x;
            shapeIntercept[row][0] = 1.0;
            response[row] = mean * noise;
        }
        DistributionalResult result = DistributionalModel.fit(response,
            List.of(PenalizedPredictor.linear(meanDesign),
                PenalizedPredictor.linear(shapeIntercept)),
            new GammaMeanShapeFamily(), DistributionalOptions.defaults(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(0.4, result.parameter("mu").coefficients()[0], 0.04);
        assertEquals(0.3, result.parameter("mu").coefficients()[1], 0.05);
        assertTrue(result.parameter("shape").fittedValues()[0] > 2.0);
    }

    @Test
    void betaMeanAndPrecisionFitJointly() {
        int observations = 200;
        double[] response = new double[observations];
        double[][] meanDesign = new double[observations][2];
        double[][] precisionIntercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            double x = -1.0 + 2.0 * row / (observations - 1.0);
            double mean = logistic(-0.2 + 0.7 * x);
            double perturbation = 0.11 * Math.sqrt(2.0)
                * Math.sin(23.0 * row + 0.7);
            meanDesign[row][0] = 1.0;
            meanDesign[row][1] = x;
            precisionIntercept[row][0] = 1.0;
            response[row] = Math.max(0.01, Math.min(0.99,
                mean + perturbation));
        }
        DistributionalResult result = DistributionalModel.fit(response,
            List.of(PenalizedPredictor.linear(meanDesign),
                PenalizedPredictor.linear(precisionIntercept)),
            new BetaMeanPrecisionFamily(), DistributionalOptions.defaults(),
            BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(-0.2, result.parameter("mu").coefficients()[0], 0.08);
        assertEquals(0.7, result.parameter("mu").coefficients()[1], 0.12);
        assertTrue(result.parameter("precision").fittedValues()[0] > 5.0);
    }

    private static double logistic(double value) {
        return 1.0 / (1.0 + Math.exp(-value));
    }
}
