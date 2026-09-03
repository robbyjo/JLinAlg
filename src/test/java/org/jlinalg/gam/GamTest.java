/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.mixed.LinearMixedModel;
import org.jlinalg.mixed.LinearMixedModelResult;
import org.jlinalg.mixed.RandomEffectTerm;
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

    @Test
    void coefficientSpaceFitMatchesDenseObservationSpaceReml() {
        int observations = 48;
        double[] x = sequence(observations, -1.0, 2.0);
        double[] response = new double[observations];
        double[] intercept = new double[observations];
        java.util.Arrays.fill(intercept, 1.0);
        for (int row = 0; row < observations; row++) {
            response[row] = 0.8 + 0.3 * x[row]
                + Math.sin(2.3 * x[row]) + 0.04 * Math.cos(11.0 * row);
        }
        PSplineTerm smooth = PSplineTerm.of("s(x)", x, 9);
        PSplineMixedModelCompiler.Compiled compiled;
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            compiled = PSplineMixedModelCompiler.compile(
                intercept, observations, 1, List.of(smooth),
                context.backend());
        }
        PSplineMixedModelCompiler.Term term = compiled.terms().get(0);
        List<String> names = java.util.stream.IntStream
            .range(0, term.randomColumns())
            .mapToObj(index -> "s(x).pen" + (index + 1)).toList();
        RandomEffectTerm random = RandomEffectTerm.of("s(x)",
            term.randomDesign(), observations, term.randomColumns(), names);
        LinearMixedModelResult dense = LinearMixedModel.fit(
            response, compiled.fixedDesign(), observations,
            compiled.fixedColumns(), List.of(random),
            RemlOptions.defaults(), BackendPolicy.CPU);
        GamResult coefficientSpace = Gam.fitGaussian(
            response, intercept, observations, 1, List.of(smooth),
            RemlOptions.defaults(), BackendPolicy.CPU);

        assertArrayEquals(dense.beta(),
            coefficientSpace.mixedModel().beta(), 2e-6);
        assertArrayEquals(dense.reml().varianceComponents(),
            coefficientSpace.mixedModel().reml().varianceComponents(), 2e-5);
        assertArrayEquals(dense.fittedValues(),
            coefficientSpace.fittedValues(), 2e-6);
        assertArrayEquals(dense.randomEffects("s(x)")
                .predictionErrorVariances(),
            coefficientSpace.mixedModel().randomEffects("s(x)")
                .predictionErrorVariances(), 2e-6);
    }

    @Test
    void explicitSatterthwaiteRequestRetainsDenseReferencePath() {
        int observations = 32;
        double[] x = sequence(observations, 0.0, 1.0);
        double[] response = new double[observations];
        double[] intercept = new double[observations];
        java.util.Arrays.fill(intercept, 1.0);
        for (int row = 0; row < observations; row++)
            response[row] = 1.0 + Math.sin(2.0 * Math.PI * x[row])
                + 0.05 * Math.cos(7.0 * row);
        RemlOptions options = RemlOptions.builder()
            .degreesOfFreedomMethod(DegreesOfFreedomMethod.SATTERTHWAITE)
            .build();

        GamResult result = Gam.fitGaussian(response, intercept,
            observations, 1, List.of(PSplineTerm.of("s(x)", x, 8)),
            options, BackendPolicy.CPU);

        assertEquals(DegreesOfFreedomMethod.SATTERTHWAITE,
            result.mixedModel().associationStatistics()
                .degreesOfFreedomMethod());
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
