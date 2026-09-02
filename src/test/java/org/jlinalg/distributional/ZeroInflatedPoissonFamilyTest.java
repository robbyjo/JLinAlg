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

class ZeroInflatedPoissonFamilyTest {
    @Test
    void separatesCountMeanFromStructuralZeroProbability() {
        int observations = 300;
        double[] response = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            intercept[row][0] = 1.0;
            if (row % 4 == 0) {
                response[row] = 0.0;
            } else {
                int pattern = (row * 17) % 11;
                response[row] = pattern < 2 ? 0.0
                    : pattern < 5 ? 1.0
                    : pattern < 8 ? 2.0
                    : pattern < 10 ? 3.0 : 4.0;
            }
        }
        PenalizedPredictor predictor = PenalizedPredictor.linear(intercept);
        DistributionalResult result = DistributionalModel.fit(response,
            List.of(predictor, predictor), new ZeroInflatedPoissonFamily(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        double countMean = result.parameter("mu").fittedValues()[0];
        double zeroProbability = result.parameter("zeroProbability")
            .fittedValues()[0];
        double observedMean = java.util.Arrays.stream(response)
            .average().orElseThrow();
        assertTrue(countMean > observedMean);
        assertTrue(zeroProbability > 0.05 && zeroProbability < 0.8);
    }
}
