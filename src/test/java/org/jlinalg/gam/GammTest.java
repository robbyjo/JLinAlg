/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
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

    private static double meanSquaredError(double[] observed, double[] fitted) {
        double sum = 0.0;
        for (int index = 0; index < observed.length; index++) {
            double residual = observed[index] - fitted[index];
            sum += residual * residual;
        }
        return sum / observed.length;
    }
}
