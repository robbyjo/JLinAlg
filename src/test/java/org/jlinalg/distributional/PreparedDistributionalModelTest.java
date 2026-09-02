/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;
import org.junit.jupiter.api.Test;

final class PreparedDistributionalModelTest {
    @Test
    void warmRefitUsesNoMoreIterationsForNearbyResponse() {
        int observations = 160;
        double[][] design = new double[observations][2];
        double[] first = new double[observations];
        double[] second = new double[observations];
        for (int row = 0; row < observations; row++) {
            double x = -1.0 + 2.0 * row / (observations - 1.0);
            design[row][0] = 1.0;
            design[row][1] = x;
            first[row] = 0.7 + 1.1 * x + 0.5 * Math.sin(17.0 * row);
            second[row] = first[row] + 0.01 * Math.cos(13.0 * row);
        }
        PenalizedPredictor predictor = PenalizedPredictor.linear(design);
        PreparedDistributionalModel prepared = new PreparedDistributionalModel(
            List.of(predictor, PenalizedPredictor.linear(intercept(observations))),
            DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
        DistributionalResult initial = prepared.fit(first);
        DistributionalResult cold = prepared.fit(second);
        DistributionalResult warm = prepared.refit(initial, second);
        assertTrue(warm.converged(), warm.convergenceMessage());
        assertTrue(warm.iterations() <= cold.iterations());
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
