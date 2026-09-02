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
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.ModelTable;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

class GeeAdvancedTest {
    @Test
    void preparedParallelAndSequentialFitsAgree() {
        Sample sample = sample(40, 4);
        double[] flat = flatten(sample.design());
        GeeOptions preparation = GeeOptions.builder()
            .correlation(GeeCorrelation.EXCHANGEABLE).build();
        PreparedGeeData prepared = Gee.prepare(sample.response(), flat,
            sample.response().length, 2, sample.cluster(), sample.wave(),
            null, null, preparation, GlmFamilies.gaussian());

        GeeResult sequential = Gee.fit(prepared, GlmFamilies.gaussian(),
            preparation, BackendPolicy.CPU);
        GeeResult parallel = Gee.fit(prepared, GlmFamilies.gaussian(),
            preparation.toBuilder().parallelism(4).parallelThreshold(1).build(),
            BackendPolicy.CPU);

        assertArrayEquals(sequential.coefficients(), parallel.coefficients(), 1e-11);
        assertArrayEquals(sequential.robustCovariance(),
            parallel.robustCovariance(), 1e-11);
        assertEquals(sample.response().length, prepared.observations());
    }

    @Test
    void projectsIndefiniteFixedCorrelationAndReportsDiagnostics() {
        Sample sample = sample(24, 3);
        double[][] indefinite = {
            {1.0, 0.60, 0.60},
            {0.60, 1.0, -0.30},
            {0.60, -0.30, 1.0}
        };
        GeeResult result = Gee.fit(sample.response(), sample.design(),
            sample.cluster(), sample.wave(), GlmFamilies.gaussian(), null, null,
            GeeOptions.builder().correlation(GeeCorrelation.FIXED)
                .fixedAssociation(indefinite).build(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertEquals(result.clusters(), result.diagnostics().clusterIds().length);
        assertEquals(result.clusters() * result.parameters(),
            result.diagnostics().oneStepDeletedCoefficients().length);
        assertTrue(Arrays.stream(result.standardizedResiduals()).allMatch(Double::isFinite));
    }

    @Test
    void smallSampleCovariancesClusterTAndExactDeletionAreAvailable() {
        Sample sample = sample(12, 3);
        GeeResult result = Gee.fit(sample.response(), sample.design(),
            sample.cluster(), sample.wave(), GlmFamilies.gaussian(), null, null,
            GeeOptions.builder().covariance(GeeCovariance.JACKKNIFE)
                .inference(GeeInference.CLUSTER_T)
                .exactClusterDeletion(true).build(), BackendPolicy.CPU);

        assertTrue(Arrays.stream(result.kauermannCarrollCovariance())
            .allMatch(Double::isFinite));
        assertTrue(Arrays.stream(result.fayGraubardCovariance())
            .allMatch(Double::isFinite));
        assertTrue(Arrays.stream(result.jackknifeCovariance())
            .allMatch(Double::isFinite));
        assertEquals(10.0, result.degreesOfFreedom());
        assertTrue(result.diagnostics().hasExactDeletionFits());
        assertEquals(GeeInference.CLUSTER_T, result.inference());
    }

    @Test
    void formulaDeclarationsPredictionAndTidyOutputWork() {
        Sample sample = sample(16, 3);
        int n = sample.response().length;
        double[] x = new double[n];
        double[] id = new double[n];
        double[] wave = new double[n];
        for (int row = 0; row < n; row++) {
            x[row] = sample.design()[row][1];
            id[row] = sample.cluster()[row];
            wave[row] = sample.wave()[row];
        }
        ModelTable table = ModelTable.builder(n)
            .numeric("y", sample.response()).numeric("x", x)
            .numeric("id", id).numeric("visit", wave).build();
        GeeResult result = GeeFormula.fit(
            "y ~ x + cluster(id) + wave(visit)", table,
            GlmFamilies.gaussian(), GeeOptions.defaults(), BackendPolicy.CPU);

        GeePrediction[] prediction = result.predict(
            new double[][] {{1.0, 1.5}}, GlmFamilies.gaussian());
        assertEquals(1, prediction.length);
        assertTrue(Double.isFinite(prediction[0].mean()));
        assertEquals(result.parameters(), result.coefficientTable().length);
        assertArrayEquals(result.responseResiduals(),
            result.residuals(GeeResidualType.RESPONSE));
        String json = result.toJson();
        assertTrue(json.startsWith("{\"family\":"));
        assertTrue(json.contains("\"coefficients\":["));
        assertTrue(json.contains("\"selectedBackend\":"));
        assertTrue(json.endsWith("}"));
    }

    @Test
    void quasiLikelihoodAndConstrainedNuisanceLinksAreFinite() {
        int clusters = 30;
        int size = 3;
        double[] response = new double[clusters * size];
        double[][] design = new double[response.length][2];
        int[] id = new int[response.length];
        int[] wave = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < size; visit++) {
                int row = cluster * size + visit;
                response[row] = (cluster + 2 * visit) % 4;
                design[row][0] = 1.0;
                design[row][1] = visit;
                id[row] = cluster;
                wave[row] = visit;
            }
        }
        GeeResult result = Gee.fit(response, design, id, wave,
            GlmFamilies.poisson(), null, null,
            GeeOptions.builder().correlation(GeeCorrelation.USER_DEFINED)
                .correlationDesign(new double[][] {{1.0}, {1.0}, {1.0}})
                .associationLink(GeeParameterLink.FISHER_Z).build(),
            BackendPolicy.CPU);
        assertTrue(Double.isFinite(result.criteria().quasiLikelihood()));
        assertTrue(Double.isFinite(result.criteria().qic()));
    }

    @Test
    void singletonClustersAndNearBoundaryCorrelationRemainStable() {
        Sample sample = sample(20, 3);
        double[][] correlation = {
            {1.0, 0.979, 0.958},
            {0.979, 1.0, 0.979},
            {0.958, 0.979, 1.0}
        };
        GeeResult result = Gee.fit(sample.response(), sample.design(),
            sample.cluster(), sample.wave(), GlmFamilies.gaussian(), null, null,
            GeeOptions.builder().correlation(GeeCorrelation.FIXED)
                .fixedAssociation(correlation).build(), BackendPolicy.CPU);
        assertTrue(Arrays.stream(result.coefficients()).allMatch(Double::isFinite));

        assertThrows(IllegalArgumentException.class, () -> Gee.fit(
            new double[] {1, 2, 3, 2, 3, 4},
            new double[][] {{1, 0}, {1, 1}, {1, 2},
                {1, 0}, {1, 1}, {1, 2}},
            new int[] {0, 0, 0, 1, 1, 1},
            new int[] {0, 1, 2, 0, 1, 2},
            GlmFamilies.gaussian(), null, null,
            GeeOptions.builder().inference(GeeInference.CLUSTER_T).build(),
            BackendPolicy.CPU));
    }

    private static Sample sample(int clusters, int size) {
        double[] response = new double[clusters * size];
        double[][] design = new double[response.length][2];
        int[] id = new int[response.length];
        int[] wave = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < size; visit++) {
                int row = cluster * size + visit;
                design[row][0] = 1.0;
                design[row][1] = visit;
                response[row] = 1.0 + 0.4 * visit
                    + (cluster % 5 - 2) * 0.08
                    + ((cluster + visit) % 3 - 1) * 0.03;
                id[row] = cluster;
                wave[row] = visit;
            }
        }
        return new Sample(response, design, id, wave);
    }

    private static double[] flatten(double[][] matrix) {
        double[] result = new double[matrix.length * matrix[0].length];
        for (int row = 0; row < matrix.length; row++) {
            System.arraycopy(matrix[row], 0, result,
                row * matrix[row].length, matrix[row].length);
        }
        return result;
    }

    private record Sample(
            double[] response, double[][] design,
            int[] cluster, int[] wave) { }
}
