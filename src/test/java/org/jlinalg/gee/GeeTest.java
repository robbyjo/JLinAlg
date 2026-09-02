/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.stream.IntStream;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.CompiledFormula;
import org.jlinalg.formula.Formula;
import org.jlinalg.formula.ModelTable;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glm.GlmResult;
import org.jlinalg.model.MissingDataPolicy;
import org.junit.jupiter.api.Test;

class GeeTest {
    @Test
    void independenceEstimateMatchesGaussianGlm() {
        double[] response = {1, 3, 2, 5, 4, 8, 6, 9};
        double[][] design = design(response.length);
        int[] cluster = {0, 0, 1, 1, 2, 2, 3, 3};

        GeeResult gee = Gee.fit(response, design, cluster, null,
            GlmFamilies.gaussian(), null, null,
            GeeOptions.defaults(), BackendPolicy.CPU);
        GlmResult glm = Glm.fit(response, design, GlmFamilies.gaussian(),
            null, null, GlmOptions.defaults(), BackendPolicy.CPU);

        assertTrue(gee.converged(), gee.convergenceMessage());
        assertArrayEquals(glm.coefficients(), gee.coefficients(), 1e-9);
        assertEquals(4, gee.clusters());
        assertEquals(2, gee.minimumClusterSize());
        assertTrue(Double.isFinite(gee.criteria().qic()));
        assertTrue(Arrays.stream(gee.robustCovariance()).allMatch(Double::isFinite));
    }

    @Test
    void builtInCorrelationStructuresFitUnequalClusters() {
        Synthetic sample = syntheticGaussian();
        for (GeeCorrelation correlation : new GeeCorrelation[] {
                GeeCorrelation.EXCHANGEABLE, GeeCorrelation.AR1,
                GeeCorrelation.M_DEPENDENT, GeeCorrelation.TOEPLITZ,
                GeeCorrelation.UNSTRUCTURED}) {
            GeeOptions options = GeeOptions.builder()
                .correlation(correlation).dependenceOrder(2).build();
            GeeResult result = Gee.fit(sample.response(), sample.design(),
                sample.cluster(), sample.repeated(), GlmFamilies.gaussian(),
                null, null, options, BackendPolicy.CPU);
            assertTrue(Arrays.stream(result.coefficients()).allMatch(Double::isFinite),
                correlation.toString());
            assertTrue(Arrays.stream(result.associationParameters())
                .allMatch(Double::isFinite), correlation.toString());
        }
    }

    @Test
    void fixedAndUserDefinedAssociationAreSupported() {
        Synthetic sample = syntheticGaussian();
        double[][] fixed = {
            {1.0, 0.2, 0.1},
            {0.2, 1.0, 0.2},
            {0.1, 0.2, 1.0}
        };
        GeeResult fixedFit = Gee.fit(sample.response(), sample.design(),
            sample.cluster(), sample.repeated(), GlmFamilies.gaussian(),
            null, null, GeeOptions.builder()
                .correlation(GeeCorrelation.FIXED)
                .fixedAssociation(fixed).build(), BackendPolicy.CPU);
        assertArrayEquals(new double[] {0.2, 0.1, 0.2},
            fixedFit.associationParameters(), 1e-12);

        double[][] pairDesign = {{1.0}, {1.0}, {1.0}};
        GeeResult userFit = Gee.fit(sample.response(), sample.design(),
            sample.cluster(), sample.repeated(), GlmFamilies.gaussian(),
            null, null, GeeOptions.builder()
                .correlation(GeeCorrelation.USER_DEFINED)
                .correlationDesign(pairDesign).build(), BackendPolicy.CPU);
        assertEquals(1, userFit.associationParameters().length);
        assertTrue(Double.isFinite(userFit.associationParameters()[0]));
    }

    @Test
    void covarianceCorrectionsAndAdjustedMethodsAreAvailable() {
        Synthetic sample = syntheticGaussian();
        for (GeeCovariance covariance : GeeCovariance.values()) {
            GeeResult result = Gee.fit(sample.response(), sample.design(),
                sample.cluster(), sample.repeated(), GlmFamilies.gaussian(),
                null, null, GeeOptions.builder().covariance(covariance).build(),
                BackendPolicy.CPU);
            assertTrue(Arrays.stream(result.covariance()).allMatch(Double::isFinite));
        }
        for (GeeMethod method : GeeMethod.values()) {
            GeeResult result = Gee.fit(sample.response(), sample.design(),
                sample.cluster(), sample.repeated(), GlmFamilies.gaussian(),
                null, null, GeeOptions.builder().method(method)
                    .maximumIterations(40).build(), BackendPolicy.CPU);
            assertTrue(Arrays.stream(result.coefficients()).allMatch(Double::isFinite),
                method.toString());
        }
    }

    @Test
    void observationLevelScaleRegressionIsEstimated() {
        Synthetic sample = syntheticGaussian();
        double[][] scaleDesign = new double[sample.response().length][2];
        for (int row = 0; row < scaleDesign.length; row++) {
            scaleDesign[row][0] = 1.0;
            scaleDesign[row][1] = sample.design()[row][1];
        }
        GeeResult result = Gee.fit(sample.response(), sample.design(),
            sample.cluster(), sample.repeated(), GlmFamilies.gaussian(),
            null, null, GeeOptions.builder().scaleDesign(scaleDesign).build(),
            BackendPolicy.CPU);

        assertEquals(2, result.scaleCoefficients().length);
        assertTrue(result.dispersion() > 0.0);
    }

    @Test
    void binaryOddsRatioAssociationFits() {
        int clusters = 20;
        int size = 3;
        double[] response = new double[clusters * size];
        double[][] design = new double[response.length][2];
        int[] id = new int[response.length];
        int[] wave = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < size; visit++) {
                int row = cluster * size + visit;
                response[row] = (cluster + visit + cluster / 3) % 3 == 0 ? 1.0 : 0.0;
                design[row][0] = 1.0;
                design[row][1] = visit;
                id[row] = cluster;
                wave[row] = visit;
            }
        }
        GeeOptions options = GeeOptions.builder()
            .association(GeeAssociation.ODDS_RATIO)
            .correlation(GeeCorrelation.EXCHANGEABLE).build();
        GeeResult result = Gee.fit(response, design, id, wave,
            GlmFamilies.binomial(), null, null, options, BackendPolicy.CPU);

        assertEquals(GeeAssociation.ODDS_RATIO, result.association());
        assertEquals(1, result.associationParameters().length);
        assertTrue(result.associationParameters()[0] > 0.0);
    }

    @Test
    void missingRowsAreOmittedBeforeClusterSorting() {
        double[] response = {1, Double.NaN, 3, 4, 5, 6};
        double[][] design = design(response.length);
        int[] cluster = {2, 0, 2, 0, 1, 1};
        int[] repeated = {1, 0, 0, 1, 1, 0};
        GeeOptions options = GeeOptions.builder()
            .missingDataPolicy(MissingDataPolicy.OMIT).build();

        GeeResult result = Gee.fit(response, design, cluster, repeated,
            GlmFamilies.gaussian(), null, null, options, BackendPolicy.CPU);

        assertEquals(5, result.observations());
        assertArrayEquals(new int[] {0, 2, 3, 4, 5}, result.retainedRows());
        assertEquals(1, result.omittedObservations());
    }

    @Test
    void compiledFormulaAdapterReusesModelMatrix() {
        ModelTable table = ModelTable.builder(8)
            .numeric("y", 1, 3, 2, 5, 4, 8, 6, 9)
            .numeric("x", 0, 1, 2, 3, 4, 5, 6, 7)
            .build();
        CompiledFormula formula = Formula.compile("y ~ x", table);
        int[] cluster = {0, 0, 1, 1, 2, 2, 3, 3};

        GeeResult result = GeeFormula.fit(formula, cluster, null,
            GlmFamilies.gaussian(), GeeOptions.defaults(), BackendPolicy.CPU);

        assertEquals(2, result.coefficients().length);
    }

    @Test
    void duplicateWavesWithinClusterAreRejected() {
        assertThrows(IllegalArgumentException.class, () -> Gee.fit(
            new double[] {1, 2, 3, 4}, design(4),
            new int[] {0, 0, 1, 1}, new int[] {0, 0, 0, 1},
            GlmFamilies.gaussian(), null, null,
            GeeOptions.defaults(), BackendPolicy.CPU));
    }

    private static double[][] design(int rows) {
        double[][] result = new double[rows][2];
        for (int row = 0; row < rows; row++) {
            result[row][0] = 1.0;
            result[row][1] = row;
        }
        return result;
    }

    private static Synthetic syntheticGaussian() {
        int clusters = 18;
        int observations = clusters * 3 - 3;
        double[] response = new double[observations];
        double[][] design = new double[observations][2];
        int[] id = new int[observations];
        int[] wave = new int[observations];
        int row = 0;
        for (int cluster = 0; cluster < clusters; cluster++) {
            int size = cluster < 3 ? 2 : 3;
            double clusterEffect = (cluster % 5 - 2) * 0.25;
            for (int visit = 0; visit < size; visit++) {
                response[row] = 1.5 + 0.7 * visit + clusterEffect
                    + ((cluster + visit) % 3 - 1) * 0.1;
                design[row][0] = 1.0;
                design[row][1] = visit;
                id[row] = cluster;
                wave[row] = visit;
                row++;
            }
        }
        assertEquals(observations, row);
        return new Synthetic(response, design, id, wave);
    }

    private record Synthetic(
            double[] response, double[][] design,
            int[] cluster, int[] repeated) {
    }
}
