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

final class AdditionalDistributionalFamiliesTest {
    @Test
    void adjacentCategoryOrdinalRecoversInterceptLogOdds() {
        double[] response = categories(40, 30, 20, 10);
        PenalizedPredictor intercept = PenalizedPredictor.linear(
            intercept(response.length));
        DistributionalResult fit = DistributionalModel.fit(response,
            List.of(intercept, intercept, intercept),
            DistributionalFamilies.ordinalAdjacentCategories(4),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertEquals(Math.log(40.0 / 30.0),
            fit.parameter("adjacentLogit0").coefficients()[0], 1e-8);
        assertEquals(Math.log(30.0 / 20.0),
            fit.parameter("adjacentLogit1").coefficients()[0], 1e-8);
        assertEquals(Math.log(20.0 / 10.0),
            fit.parameter("adjacentLogit2").coefficients()[0], 1e-8);
    }

    @Test
    void hurdlePoissonSeparatesZerosFromPositiveCounts() {
        double[] response = new double[100];
        for (int row = 40; row < 70; row++) response[row] = 1.0;
        for (int row = 70; row < 90; row++) response[row] = 2.0;
        for (int row = 90; row < 100; row++) response[row] = 3.0;
        PenalizedPredictor intercept = PenalizedPredictor.linear(
            intercept(response.length));
        DistributionalResult fit = DistributionalModel.fit(response,
            List.of(intercept, intercept),
            DistributionalFamilies.hurdlePoisson(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
        assertTrue(fit.converged(), fit.convergenceMessage());
        assertEquals(0.4, fit.parameter("zeroProbability").fittedValues()[0],
            1e-8);
        assertTrue(fit.parameter("mu").fittedValues()[0] > 1.0);
        assertTrue(fit.parameter("mu").fittedValues()[0] < 1.3);
    }

    @Test
    void negativeBinomialScoresAndInformationAreFinite() {
        DistributionalFamily family =
            DistributionalFamilies.negativeBinomialMeanDispersion();
        double[] score = new double[2];
        double[] information = new double[4];
        family.derivatives(4.0, new double[] {2.5, 1.7}, score, information);
        for (double value : score) assertTrue(Double.isFinite(value));
        for (double value : information) assertTrue(Double.isFinite(value));
        assertTrue(information[0] > 0.0);
        assertTrue(information[3] > 0.0);
    }

    private static double[] categories(int... counts) {
        int total = 0;
        for (int count : counts) total += count;
        double[] result = new double[total];
        int row = 0;
        for (int category = 0; category < counts.length; category++) {
            for (int count = 0; count < counts[category]; count++) {
                result[row++] = category;
            }
        }
        return result;
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
