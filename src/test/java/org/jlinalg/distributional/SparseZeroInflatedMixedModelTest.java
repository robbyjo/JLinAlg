/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
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
    void fullObservedHessianMatchesTwoDimensionalFiniteDifferences() {
        double step = 1e-4;
        for (boolean negativeBinomial : new boolean[] {false, true}) {
            for (double response : new double[] {0.0, 3.0}) {
                double count = 0.7;
                double zero = -0.8;
                double[] center = SparseZeroInflatedMixedModel
                    .likelihoodDerivatives(response, count, zero, 2.3,
                        negativeBinomial);
                double[] countPlus = SparseZeroInflatedMixedModel
                    .likelihoodDerivatives(response, count + step, zero, 2.3,
                        negativeBinomial);
                double[] countMinus = SparseZeroInflatedMixedModel
                    .likelihoodDerivatives(response, count - step, zero, 2.3,
                        negativeBinomial);
                double[] zeroPlus = SparseZeroInflatedMixedModel
                    .likelihoodDerivatives(response, count, zero + step, 2.3,
                        negativeBinomial);
                double[] zeroMinus = SparseZeroInflatedMixedModel
                    .likelihoodDerivatives(response, count, zero - step, 2.3,
                        negativeBinomial);
                assertEquals(-(zeroPlus[0] - 2.0 * center[0] + zeroMinus[0])
                    / (step * step), center[4], 3e-6);
                assertEquals(-(countPlus[2] - countMinus[2]) / (2.0 * step),
                    center[5], 3e-7);
                assertEquals(-(zeroPlus[1] - zeroMinus[1]) / (2.0 * step),
                    center[5], 3e-7);
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
    void selectableOuterOptimizersAgreeOnGroupedZip() {
        Fixture fixture = fixture(false);
        ZeroInflatedMixedOptions base = controls();
        ZeroInflatedMixedResult bfgs = SparseZeroInflatedMixedModel.fitPoisson(
            fixture.response(), fixture.countFixed(), 2,
            fixture.intercept(), 1, List.of(fixture.groupTerm()), null, null,
            withOptimizer(base, ZeroInflatedOuterOptimizer.BOUNDED_BFGS),
            BackendPolicy.CPU);
        ZeroInflatedMixedResult bobyqa = SparseZeroInflatedMixedModel.fitPoisson(
            fixture.response(), fixture.countFixed(), 2,
            fixture.intercept(), 1, List.of(fixture.groupTerm()), null, null,
            withOptimizer(base, ZeroInflatedOuterOptimizer.BOBYQA),
            BackendPolicy.CPU);

        assertTrue(bfgs.converged(), bfgs.convergenceMessage());
        assertTrue(bobyqa.converged(), bobyqa.convergenceMessage());
        assertEquals(bobyqa.marginalLogLikelihood(),
            bfgs.marginalLogLikelihood(), 1e-4);
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
    void outerCachesRefreshForDispersionOffsetsAndParallelGradients() {
        Fixture fixture = fixture(true);
        int rows = fixture.response().length;
        double[] count = new double[rows * 3];
        double[] offsets = new double[rows];
        for (int row = 0; row < rows; row++) {
            count[3 * row] = 1.0;
            count[3 * row + 1] = fixture.countFixed()[2 * row + 1];
            count[3 * row + 2] = Math.sin(row * 0.73);
            offsets[row] = 0.12 * Math.cos(row * 0.31);
        }
        // Eight/nine outer parameters exercise the actual parallel-gradient path.
        for (boolean zeroRandom : new boolean[] {false, true}) {
            List<RandomEffectTerm> zeroTerms = zeroRandom
                ? List.of(fixture.groupTerm()) : List.of();
            try (SparseZeroInflatedMixedModel.Prepared serial =
                    SparseZeroInflatedMixedModel.prepareNegativeBinomial(rows,
                        List.of(fixture.groupTerm()), null, zeroTerms, null,
                        List.of(), cacheTestControls(1), BackendPolicy.CPU);
                 SparseZeroInflatedMixedModel.Prepared parallel =
                    SparseZeroInflatedMixedModel.prepareNegativeBinomial(rows,
                        List.of(fixture.groupTerm()), null, zeroTerms, null,
                        List.of(), cacheTestControls(4), BackendPolicy.CPU)) {
                for (int pass = 0; pass < 2; pass++) {
                    double[] response = fixture.response().clone();
                    double[] dispersion = fixture.countFixed().clone();
                    if (pass == 1) {
                        // Reuse symbolic factors with different data/design/offsets.
                        for (int row = 0; row < rows; row++) {
                            if (row % 7 == 0) response[row] += 1.0;
                            dispersion[2 * row + 1] = Math.cos(row * 0.17);
                            offsets[row] += 0.02;
                        }
                    }
                    ZeroInflatedMixedResult expected = serial.fit(response,
                        count, 3, fixture.countFixed(), 2, dispersion, 2, offsets);
                    ZeroInflatedMixedResult actual = parallel.fit(response,
                        count, 3, fixture.countFixed(), 2, dispersion, 2, offsets);
                    System.out.printf("ZI cache zeroRandom=%s pass=%d serial=%s calls=%d modes=%d LL=%.12f parallel=%s calls=%d%n",
                        zeroRandom,pass,expected.converged(),expected.objectiveEvaluations(),expected.modeIterations(),expected.marginalLogLikelihood(),actual.converged(),actual.objectiveEvaluations());
                    assertTrue(expected.converged(), expected.convergenceMessage());
                    assertTrue(actual.converged(), actual.convergenceMessage());
                    assertTrue(expected.optimizationStarts() >= 5);
                    assertEquals(expected.globalOptimumCertified(),
                        actual.globalOptimumCertified());
                    assertTrue(expected.globalOptimumCertified(),
                        "selected likelihood must be best across multiple stationary starts");
                    assertArrayEquals(expected.outerParameterEstimates(),
                        actual.outerParameterEstimates(), 1e-8);
                    assertArrayEquals(expected.fittedMeans(), actual.fittedMeans(), 1e-8);
                    assertEquals(expected.objectiveEvaluations(), actual.objectiveEvaluations());
                    assertEquals(expected.marginalLogLikelihood(),
                        actual.marginalLogLikelihood(), 1e-7);
                    assertUncachedGroupedLikelihood(actual, response, count,
                        fixture.countFixed(), dispersion, offsets, zeroRandom);
                }
            }
        }
    }

    private static ZeroInflatedMixedOptions cacheTestControls(int threads) {
        return new ZeroInflatedMixedOptions(2000, 100, 1e-6, 0.4,
            1e-6, 100.0, 1e-4, 1e4, 20.0, null,
            ZeroInflatedOuterOptimizer.BOUNDED_BFGS, threads);
    }

    @Test void profileCannotCertifyAnUnconvergedRandomMode() {
        Fixture f = fixture(false);
        var fit = SparseZeroInflatedMixedModel.fitPoisson(f.response(), f.countFixed(), 2,
            f.intercept(), 1, List.of(f.groupTerm()), null, null, controls(), BackendPolicy.CPU);
        var limited = new ZeroInflatedMixedOptions(500, 1, 1e-6, .4, 1e-6, 100, 1e-4, 1e4, 20, null);
        try (var prepared = SparseZeroInflatedMixedModel.preparePoisson(f.response().length,
                List.of(f.groupTerm()), null, List.of(), null, List.of(), limited, BackendPolicy.CPU)) {
            var profile = prepared.profile(fit, f.response(), f.countFixed(), 2, f.intercept(), 1,
                null, 0, null, 0, new double[] {fit.countCoefficients()[0] + .5});
            assertFalse(profile.converged()[0]);
            assertTrue(Double.isNaN(profile.logLikelihoods()[0]));
        }
    }

    /** Independently reconstruct a block-diagonal Laplace likelihood, uncached. */
    private static void assertUncachedGroupedLikelihood(
            ZeroInflatedMixedResult fit, double[] response, double[] count,
            double[] zero, double[] dispersion, double[] offsets,
            boolean zeroRandom) {
        double[] beta = fit.countCoefficients();
        double[] gamma = fit.zeroCoefficients();
        double[] delta = fit.dispersionCoefficients();
        double[] u = fit.randomEffects("group");
        double[] v = zeroRandom ? fit.randomEffects("zero:group") : new double[u.length];
        double[] variances = fit.varianceComponents();
        double[] fittedSizes = fit.sizes();
        double[] fittedCounts = fit.conditionalCountMeans();
        double[] fittedZeros = fit.structuralZeroProbabilities();
        double likelihood = 0.0;
        for (int group = 0; group < u.length; group++) {
            double hcc = 1.0 / variances[0];
            double hzz = zeroRandom ? 1.0 / variances[1] : 1.0;
            double hcz = 0.0;
            for (int within = 0; within < 20; within++) {
                int row = 20 * group + within;
                double fixedCount = 0.0;
                for (int column = 0; column < 3; column++)
                    fixedCount += count[3 * row + column] * beta[column];
                double eta = offsets[row] + fixedCount + u[group];
                double zeta = zero[2 * row] * gamma[0]
                    + zero[2 * row + 1] * gamma[1] + v[group];
                double size = Math.exp(dispersion[2 * row] * delta[0]
                    + dispersion[2 * row + 1] * delta[1]);
                assertEquals(size, fittedSizes[row], 1e-12);
                assertEquals(Math.exp(eta), fittedCounts[row], 1e-12);
                assertEquals(1.0 / (1.0 + Math.exp(-zeta)), fittedZeros[row], 1e-12);
                double[] uncached = SparseZeroInflatedMixedModel.likelihoodDerivatives(
                    response[row], eta, zeta, size, true);
                likelihood += uncached[0];
                hcc += uncached[3];
                if (zeroRandom) {
                    hzz += uncached[4];
                    hcz += uncached[5];
                }
            }
            likelihood -= 0.5 * (u[group] * u[group] / variances[0]
                + Math.log(variances[0]) + Math.log(hcc * hzz - hcz * hcz));
            if (zeroRandom)
                likelihood -= 0.5 * (v[group] * v[group] / variances[1]
                    + Math.log(variances[1]));
        }
        assertEquals(likelihood, fit.marginalLogLikelihood(), 1e-9);
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

    @Test
    void preparedFitSupportsIndependentZeroProcessRandomEffects() {
        Fixture fixture = fixture(false);
        try (SparseZeroInflatedMixedModel.Prepared prepared =
                SparseZeroInflatedMixedModel.preparePoisson(
                    fixture.response().length,
                    List.of(fixture.groupTerm()), null,
                    List.of(fixture.groupTerm()), null, List.of(),
                    controls(), BackendPolicy.CPU)) {
            ZeroInflatedMixedResult fit = prepared.fitWithInference(
                fixture.response(),
                fixture.countFixed(), 2, fixture.intercept(), 1,
                null, 0, null);
            assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
            assertEquals(24, fit.randomCoefficientCount());
            assertEquals(2, fit.varianceComponents().length);
            assertEquals(12, fit.randomEffects("zero:group").length);
            assertEquals(prepared.sparseEquationNonzeroCount(),
                fit.sparseEquationNonzeroCount());
            assertTrue(fit.inferenceAvailable(), () ->
                "converged=" + fit.converged()
                    + " certified=" + fit.globalOptimumCertified()
                    + " starts=" + fit.optimizationStarts()
                    + " parameters="
                    + java.util.Arrays.toString(
                        fit.outerParameterEstimates())
                    + " warnings=" + fit.warnings());
            assertEquals(3, fit.standardErrors().length);
            for (double value : fit.standardErrors())
                assertTrue(Double.isFinite(value) && value >= 0.0);
            double estimate = fit.outerParameterEstimates()[0];
            ZeroInflatedProfile profile = prepared.profile(fit,
                fixture.response(), fixture.countFixed(), 2,
                fixture.intercept(), 1, null, 0, null, 0,
                new double[] {estimate});
            assertEquals(1, profile.logLikelihoods().length);
            assertTrue(Double.isFinite(profile.logLikelihoods()[0]));
            ZeroInflatedBootstrapResult bootstrap =
                prepared.parametricBootstrap(fit, fixture.countFixed(), 2,
                    fixture.intercept(), 1, null, 0, null, 1, 91L);
            assertEquals(1, bootstrap.replicates());
            assertEquals(fit.outerParameterNames(), bootstrap.parameterNames());
        }
    }

    @Test
    void correlatedPedigreeUsesGuardedTwoProcessCovariance() {
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
            response[row] = row % (id.equals("a") ? 3 : 7) == 0
                ? 0.0 : id.equals("a") ? 2.0 : 5.0;
        }
        PedigreeRandomEffectTerm pedigreeTerm =
            PedigreeRandomEffectTerm.of("pedigree", ids, pedigree);
        CorrelatedZeroInflatedRandomEffect correlated =
            CorrelatedZeroInflatedRandomEffect.pedigree(
                "additive", pedigreeTerm);
        ZeroInflatedMixedResult fit = SparseZeroInflatedMixedModel.fitPoisson(
            response, intercept, 1, intercept, 1,
            List.of(), List.of(), List.of(), List.of(), List.of(correlated),
            null, controls(), BackendPolicy.CPU);
        assertTrue(Double.isFinite(fit.marginalLogLikelihood()));
        assertEquals(2 * pedigree.size(), fit.randomCoefficientCount());
        assertEquals(2, fit.varianceComponents().length);
        assertTrue(Math.abs(fit.correlations().get("additive")) < 0.99);
        assertEquals(pedigree.size(),
            fit.randomEffects("additive:zero").length);
    }

    @Test
    void simulationRecoversZeroProcessGroupOrdering() {
        int groups = 16;
        int perGroup = 45;
        int rows = groups * perGroup;
        double[] response = new double[rows];
        double[] intercept = new double[rows];
        List<String> labels = new ArrayList<>(rows);
        Random random = new Random(8_190_517L);
        for (int row = 0; row < rows; row++) {
            int group = row / perGroup;
            intercept[row] = 1.0;
            labels.add("g" + group);
            double zeroEffect = group % 2 == 0 ? -1.15 : 1.15;
            double pi = 1.0 / (1.0 + Math.exp(-(-1.0 + zeroEffect)));
            response[row] = random.nextDouble() < pi
                ? 0.0 : poisson(random, 2.8);
        }
        RandomEffectTerm zero = RandomEffectTerm.randomIntercept(
            "group", labels);
        ZeroInflatedMixedResult fit = SparseZeroInflatedMixedModel.fitPoisson(
            response, intercept, 1, intercept, 1,
            List.of(), List.of(), List.of(zero), null, List.of(), null,
            controls(), BackendPolicy.CPU);
        double[] modes = fit.randomEffects("zero:group");
        double low = 0.0;
        double high = 0.0;
        for (int group = 0; group < groups; group++) {
            if (group % 2 == 0) low += modes[group];
            else high += modes[group];
        }
        assertTrue(high / (groups / 2.0) > low / (groups / 2.0) + 0.8);
        assertTrue(fit.varianceComponents()[0] > 0.15);
    }

    @Test
    void preparedStructureCreatesOneNumericFactorPerWorker() throws Exception {
        Fixture fixture = fixture(false);
        try (SparseZeroInflatedMixedModel.Prepared prepared =
                SparseZeroInflatedMixedModel.preparePoisson(
                    fixture.response().length, List.of(fixture.groupTerm()),
                    null, List.of(), null, List.of(), controls(),
                    BackendPolicy.CPU)) {
            ExecutorService workers = Executors.newFixedThreadPool(2);
            try {
                java.util.concurrent.Callable<ZeroInflatedMixedResult> fit = () ->
                    prepared.fit(fixture.response(), fixture.countFixed(), 2,
                        fixture.intercept(), 1, null, 0, null);
                Future<ZeroInflatedMixedResult> first = workers.submit(fit);
                Future<ZeroInflatedMixedResult> second = workers.submit(fit);
                assertTrue(Double.isFinite(first.get().marginalLogLikelihood()));
                assertTrue(Double.isFinite(second.get().marginalLogLikelihood()));
                assertEquals(2, prepared.numericFactorCount());
            } finally {
                workers.shutdownNow();
            }
        }
    }

    private static ZeroInflatedMixedOptions controls() {
        return new ZeroInflatedMixedOptions(
            350, 100, 1e-6, 0.4,
            1e-6, 100.0, 1e-4, 1e4, 20.0, null);
    }

    private static ZeroInflatedMixedOptions withOptimizer(
            ZeroInflatedMixedOptions options,
            ZeroInflatedOuterOptimizer optimizer) {
        return new ZeroInflatedMixedOptions(
            options.maximumOuterEvaluations(), options.maximumModeIterations(),
            options.relativeTolerance(), options.initialTrustRadius(),
            options.minimumVariance(), options.maximumVariance(),
            options.minimumSize(), options.maximumSize(),
            options.maximumAbsoluteCoefficient(), options.initialVariances(),
            optimizer);
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
