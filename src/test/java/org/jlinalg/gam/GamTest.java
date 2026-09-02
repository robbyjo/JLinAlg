/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class GamTest {
    @Test
    void bsplineBasisIsAPartitionOfUnity() {
        double[] x = sequence(31, 0.0, 1.0);
        PSplineTerm term = PSplineTerm.of("x", x, 9);
        double[] basis = term.design();
        for (int row = 0; row < x.length; row++) {
            double sum = 0.0;
            for (int column = 0; column < term.basisDimension(); column++) {
                double value = basis[row * term.basisDimension() + column];
                assertTrue(value >= -1e-15);
                sum += value;
            }
            assertEquals(1.0, sum, 1e-12);
        }
    }

    @Test
    void gaussianRemlRecoversANonlinearSignalAndPredictsTrainingData() {
        int observations = 80;
        double[] x = sequence(observations, 0.0, 1.0);
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            intercept[row][0] = 1.0;
            response[row] = 1.5 + Math.sin(2.0 * Math.PI * x[row])
                + 0.08 * Math.sin(17.0 * row);
        }

        GamResult result = Gam.fitGaussian(
            response, intercept, List.of(PSplineTerm.of("s(x)", x, 10)),
            RemlOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.mixedModel().reml().converged(),
            result.mixedModel().reml().convergenceMessage());
        assertTrue(result.smoothTerms().get(0).smoothingParameter() > 0.0);
        assertTrue(result.smoothTerms().get(0).effectiveDegreesOfFreedom() > 2.0);
        assertTrue(meanSquaredError(response, result.fittedValues()) < 0.03);

        double[] rowMajorIntercept = new double[observations];
        java.util.Arrays.fill(rowMajorIntercept, 1.0);
        double[] predicted = result.predict(
            rowMajorIntercept, observations, List.of(x));
        for (int row = 0; row < observations; row++) {
            assertEquals(result.fittedValues()[row], predicted[row], 1e-8);
        }
    }

    @Test
    void secondDifferencePenaltyLeavesLinearTrendUnpenalized() {
        double[] x = sequence(40, -2.0, 2.0);
        double[] response = new double[x.length];
        double[][] intercept = new double[x.length][1];
        for (int row = 0; row < x.length; row++) {
            intercept[row][0] = 1.0;
            response[row] = 3.0 - 1.25 * x[row]
                + 0.01 * Math.cos(5.0 * row);
        }
        GamResult result = Gam.fitGaussian(
            response, intercept, List.of(PSplineTerm.of("s(x)", x, 8)),
            RemlOptions.defaults(), BackendPolicy.CPU);

        assertEquals(3.0, result.parametricCoefficients()[0], 0.01);
        assertTrue(meanSquaredError(response, result.fittedValues()) < 0.001);
        double meanSmooth = java.util.Arrays.stream(
            result.smoothTerms().get(0).fittedValues()).average().orElseThrow();
        assertEquals(0.0, meanSmooth, 1e-10);
    }

    private static double[] sequence(int length, double lower, double upper) {
        double[] result = new double[length];
        for (int index = 0; index < length; index++) {
            result[index] = lower + (upper - lower) * index / (length - 1.0);
        }
        return result;
    }

    private static double meanSquaredError(double[] observed, double[] fitted) {
        double result = 0.0;
        for (int index = 0; index < observed.length; index++) {
            double residual = observed[index] - fitted[index];
            result += residual * residual;
        }
        return result / observed.length;
    }
}
