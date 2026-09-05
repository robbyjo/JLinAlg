/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.junit.jupiter.api.Test;

/** Tests the frequentist sparse Laplace zero-inflated mixed-model path. */
final class SparseZeroInflatedMixedModelTest {
    @Test
    void analyticCountDerivativesMatchFiniteDifferences() {
        double step = 1e-4;
        for (boolean negativeBinomial : new boolean[] {false, true}) {
            for (double response : new double[] {0.0, 3.0}) {
                double eta = 0.7;
                double[] center = derivatives(
                    response, eta, negativeBinomial);
                double plus = derivatives(
                    response, eta + step, negativeBinomial)[0];
                double minus = derivatives(
                    response, eta - step, negativeBinomial)[0];
                double numericalScore = (plus - minus) / (2.0 * step);
                double numericalCurvature =
                    -(plus - 2.0 * center[0] + minus) / (step * step);
                assertEquals(numericalScore, center[1], 2e-7);
                assertEquals(numericalCurvature, center[2], 2e-6);
            }
        }
    }

    @Test
    void groupedZipSeparatesStructuralZerosAndCountHeterogeneity() {
        Fixture fixture = fixture(false);
        ZeroInflatedMixedResult fit =
            SparseZeroInflatedMixedModel.fitPoisson(
                fixture.response(), fixture.countFixed(), 2,
                fixture.intercept(), 1,
                List.of(fixture.groupTerm()), null, null,
                controls(), BackendPolicy.CPU);

        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertEquals(2, fit.countCoefficients().length);
        assertEquals(1, fit.zeroCoefficients().length);
        assertTrue(fit.varianceComponents()[0] > 0.05);
        assertTrue(fit.structuralZeroProbabilities()[0] > 0.02);
        assertEquals(fixture.response().length, fit.fittedMeans().length);
        assertTrue(fit.fittedZeroProbabilities()[0]
            > fit.structuralZeroProbabilities()[0]);
    }

    @Test
    void groupedZinbEstimatesFiniteSizeAndUnconditionalPredictions() {
        Fixture fixture = fixture(true);
        ZeroInflatedMixedResult fit =
            SparseZeroInflatedMixedModel.fitNegativeBinomial(
                fixture.response(), fixture.countFixed(), 2,
                fixture.intercept(), 1, fixture.intercept(), 1,
                List.of(fixture.groupTerm()), null, null,
                controls(), BackendPolicy.CPU);

        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertEquals(1, fit.dispersionCoefficients().length);
        assertEquals(fixture.response().length, fit.sizes().length);
        assertTrue(fit.sizes()[0] > 0.1 && fit.sizes()[0] < 1000.0);
        assertTrue(fit.fittedMeans()[0] < fit.conditionalCountMeans()[0]);
    }

    @Test
    void pedigreePrecisionRetainsUnobservedAncestors() {
        Pedigree pedigree = Pedigree.of(List.of(
            PedigreeIndividual.founder("sire"),
            PedigreeIndividual.founder("dam"),
            new PedigreeIndividual("a", "sire", "dam"),
            new PedigreeIndividual("b", "sire", "dam")));
        int rows = 80;
        double[] response = new double[rows];
        double[] intercept = new double[rows];
        List<String> ids = new ArrayList<>(rows);
        for (int row = 0; row < rows; row++) {
            intercept[row] = 1.0;
            String id = row % 2 == 0 ? "a" : "b";
            ids.add(id);
            response[row] = row % 5 == 0 ? 0.0
                : id.equals("a") ? 1.0 + row % 2 : 3.0 + row % 3;
        }
        PedigreeRandomEffectTerm additive = PedigreeRandomEffectTerm.of(
            "additive", ids, pedigree);
        ZeroInflatedMixedResult fit =
            SparseZeroInflatedMixedModel.fitPoisson(
                response, intercept, 1, intercept, 1,
                List.of(additive.randomEffect()), List.of(additive.precision()),
                null, controls(), BackendPolicy.CPU);

        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
        assertEquals(pedigree.size(), fit.randomCoefficientCount());
        assertEquals(pedigree.size(), fit.randomEffects("additive").length);
        for (double value : fit.randomEffects("additive"))
            assertTrue(Double.isFinite(value));
    }

    private static ZeroInflatedMixedOptions controls() {
        return new ZeroInflatedMixedOptions(
            350, 100, 1e-6, 0.4,
            1e-6, 100.0, 1e-4, 1e4, 20.0, null);
    }

    private static double[] derivatives(
            double response, double countPredictor,
            boolean negativeBinomial) {
        return SparseZeroInflatedMixedModel.countLikelihoodDerivatives(
            response, countPredictor, -0.8, 2.3, negativeBinomial);
    }

    private static Fixture fixture(boolean negativeBinomial) {
        int groups = 12;
        int perGroup = 20;
        int rows = groups * perGroup;
        double[] response = new double[rows];
        double[] countFixed = new double[rows * 2];
        double[] intercept = new double[rows];
        List<String> labels = new ArrayList<>(rows);
        Random random = new Random(7_413_921L + (negativeBinomial ? 1 : 0));
        for (int row = 0; row < rows; row++) {
            int group = row / perGroup;
            double x = -1.0 + 2.0 * (row % perGroup) / (perGroup - 1.0);
            double randomEffect = 0.75 * Math.sin(1.1 * group);
            double mean = Math.exp(0.45 + 0.5 * x + randomEffect);
            countFixed[row * 2] = 1.0;
            countFixed[row * 2 + 1] = x;
            intercept[row] = 1.0;
            labels.add(Integer.toString(group));
            if (random.nextDouble() < 0.22) {
                response[row] = 0.0;
            } else {
                response[row] = negativeBinomial
                    ? negativeBinomial(random, mean, 2)
                    : poisson(random, mean);
            }
        }
        return new Fixture(response, countFixed, intercept,
            RandomEffectTerm.randomIntercept("group", labels));
    }

    private static int poisson(Random random, double mean) {
        double threshold = Math.exp(-mean);
        int count = 0;
        double product = 1.0;
        do {
            count++;
            product *= random.nextDouble();
        } while (product > threshold);
        return count - 1;
    }

    private static int negativeBinomial(
            Random random, double mean, int size) {
        double successProbability = size / (size + mean);
        int successes = 0;
        int failures = 0;
        while (successes < size) {
            if (random.nextDouble() < successProbability) successes++;
            else failures++;
        }
        return failures;
    }

    private record Fixture(
            double[] response, double[] countFixed, double[] intercept,
            RandomEffectTerm groupTerm) { }
}
