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

final class DistributionalDiagnosticsTest {
    @Test
    void gaussianResidualsCentilesPredictionAndComparisonAreAvailable() {
        int observations = 100;
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        double[][] linear = new double[observations][2];
        for (int row = 0; row < observations; row++) {
            double x = -1.0 + 2.0 * row / (observations - 1.0);
            response[row] = 1.0 + 0.8 * x + 0.4 * Math.sin(17.0 * row);
            intercept[row][0] = 1.0;
            linear[row][0] = 1.0;
            linear[row][1] = x;
        }
        PenalizedPredictor constant = PenalizedPredictor.linear(intercept);
        PenalizedPredictor trend = PenalizedPredictor.linear(linear);
        DistributionalResult small = DistributionalModel.fit(response,
            List.of(constant, constant),
            DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
        DistributionalResult large = DistributionalModel.fit(response,
            List.of(trend, constant),
            DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
        double[] residuals = DistributionalDiagnostics.quantileResiduals(
            response, large);
        assertEquals(observations, residuals.length);
        double[] median = DistributionalDiagnostics.gaussianCentile(large, 0.5);
        assertEquals(large.parameter("mu").fittedValues()[20], median[20], 1e-12);
        double[][] predicted = DistributionalPrediction.parameters(large,
            List.of(trend, constant),
            DistributionalFamilies.gaussianLocationScale());
        assertEquals(large.parameter("mu").fittedValues()[40], predicted[0][40], 1e-12);
        DistributionalModelComparison comparison =
            DistributionalModelComparison.compare(small, large);
        assertTrue(comparison.likelihoodRatio() > 20.0);
        assertTrue(comparison.pValue() < 1e-4);
    }
}
