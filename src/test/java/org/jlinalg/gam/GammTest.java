/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class GammTest {
    @Test
    void combinesSmoothAndGroupedCovarianceComponents() {
        int groups = 12;
        int perGroup = 5;
        int observations = groups * perGroup;
        double[] x = new double[observations];
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        List<String> labels = new ArrayList<>(observations);
        for (int row = 0; row < observations; row++) {
            int group = row / perGroup;
            x[row] = row / (observations - 1.0);
            intercept[row][0] = 1.0;
            labels.add("g" + group);
            double groupEffect = (group % 4 - 1.5) * 0.35;
            response[row] = 2.0 + Math.sin(2.0 * Math.PI * x[row])
                + groupEffect + 0.04 * Math.cos(13.0 * row);
        }

        GammResult result = Gamm.fitGaussian(response, intercept,
            List.of(PSplineTerm.of("s(x)", x, 9)),
            List.of(VarianceComponent.randomIntercept("group", labels)),
            null, RemlOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.reml().converged(),
            result.reml().convergenceMessage());
        int groupIndex = result.reml().componentNames().indexOf("group");
        assertTrue(groupIndex >= 0);
        assertTrue(result.reml().varianceComponents()[groupIndex] > 0.01);
        assertTrue(meanSquaredError(response, result.fittedValues()) < 0.02);
        assertTrue(result.randomContribution("group").length == observations);
    }

    @Test
    void preparedBatchScanMatchesDenseObservationSpaceGamm() {
        int groups = 8;
        int perGroup = 4;
        int observations = groups * perGroup;
        double[] x = new double[observations];
        double[] response = new double[observations];
        double[] intercept = new double[observations];
        double[][] interceptMatrix = new double[observations][1];
        List<String> labels = new ArrayList<>(observations);
        for (int row = 0; row < observations; row++) {
            int group = row / perGroup;
            x[row] = row / (observations - 1.0);
            intercept[row] = 1.0;
            interceptMatrix[row][0] = 1.0;
            labels.add("g" + group);
            response[row] = 1.5 + Math.sin(2.0 * Math.PI * x[row])
                + 0.18 * (group % 3 - 1)
                + 0.03 * Math.cos(7.0 * row);
        }
        RemlOptions options = RemlOptions.defaults();
        GammResult dense = Gamm.fitGaussian(response, interceptMatrix,
            List.of(PSplineTerm.of("s(x)", x, 8)),
            List.of(VarianceComponent.randomIntercept("batch", labels)),
            null, options, BackendPolicy.CPU);
        RandomEffectTerm batch =
            RandomEffectTerm.randomIntercept("batch", labels);
        try (PreparedGammPredictorScan scan =
                new PreparedGammPredictorScan(response, intercept,
                    observations, 1, List.of(batch),
                    List.of(SparsePrecisionMatrix.identity(
                        batch.coefficients())),
                    8, options, BackendPolicy.CPU)) {
            GammScanResult sparse = scan.fit("s(x)", x);
            assertEquals(dense.reml().logLikelihood(),
                sparse.reml().logLikelihood(), 2e-5);
            assertEquals(dense.smoothTerms().get(0)
                    .effectiveDegreesOfFreedom(),
                sparse.smoothTerm().effectiveDegreesOfFreedom(), 2e-4);
            assertEquals(meanSquaredError(
                    dense.fittedValues(), sparse.fittedValues()),
                0.0, 2e-8);
            assertTrue(sparse.randomContribution("batch").length
                == observations);
        }
    }

    private static double meanSquaredError(double[] observed, double[] fitted) {
        double sum = 0.0;
        for (int index = 0; index < observed.length; index++) {
            double residual = observed[index] - fitted[index];
            sum += residual * residual;
        }
        return sum / observed.length;
    }
}
