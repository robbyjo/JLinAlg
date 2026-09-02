/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;
import org.jlinalg.gam.QuadraticPenalizedPredictor;
import org.jlinalg.gam.QuadraticSmoothTerm;
import org.jlinalg.gam.TensorProductPSplineTerm;
import org.junit.jupiter.api.Test;

final class TensorDistributionalModelTest {
    @Test
    void tensorSmoothCapturesBivariateMeanSurface() {
        int observations = 144;
        double[] first = new double[observations];
        double[] second = new double[observations];
        double[] response = new double[observations];
        for (int row = 0; row < observations; row++) {
            first[row] = (row % 12) / 11.0;
            second[row] = (row / 12) / 11.0;
            response[row] = 1.2 + Math.sin(2.0 * Math.PI * first[row])
                + 0.65 * Math.cos(2.0 * Math.PI * second[row])
                + 0.04 * Math.sin(13.0 * row);
        }
        QuadraticSmoothTerm surface = TensorProductPSplineTerm.of(
            "te(x,z)", first, second, 7, 7);
        PenalizedPredictor mean = QuadraticPenalizedPredictor.compile(
            intercept(observations), List.of(surface),
            List.of(new double[] {0.2, 0.2}), BackendPolicy.CPU);
        DistributionalResult fit = DistributionalModel.fit(response,
            List.of(mean, PenalizedPredictor.linear(intercept(observations))),
            DistributionalFamilies.gaussianLocationScale(),
            DistributionalOptions.defaults(), BackendPolicy.CPU);
        assertTrue(fit.converged(), fit.convergenceMessage());
        double[] fitted = fit.parameter("mu").fittedValues();
        double rss = 0.0;
        double nullRss = 0.0;
        for (int row = 0; row < observations; row++) {
            double error = response[row] - fitted[row];
            rss += error * error;
            double centered = response[row] - 1.2;
            nullRss += centered * centered;
        }
        assertTrue(rss < 0.08 * nullRss, "tensor surface should explain the signal");
    }

    private static double[][] intercept(int observations) {
        double[][] result = new double[observations][1];
        for (double[] row : result) row[0] = 1.0;
        return result;
    }
}
