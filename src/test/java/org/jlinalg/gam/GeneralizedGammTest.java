/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class GeneralizedGammTest {
    @Test
    void poissonModelCombinesSmoothAndGroupedCovariance() {
        int groups = 10;
        int perGroup = 8;
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
            double eta = 1.1 + 0.6 * Math.sin(2.0 * Math.PI * x[row])
                + (group % 3 - 1.0) * 0.2;
            double fraction = ((row * 29) % 97) / 97.0;
            response[row] = Math.floor(Math.exp(eta) + fraction);
        }

        GeneralizedGammResult result = GeneralizedGamm.fit(
            response, intercept, List.of(PSplineTerm.of("s(x)", x, 8)),
            GlmFamilies.poisson(),
            List.of(VarianceComponent.randomIntercept("group", labels)),
            null, null, GlmmPqlOptions.defaults(), BackendPolicy.CPU);

        assertTrue(result.converged(), result.convergenceMessage());
        assertTrue(result.workingModel().conditionalDeviance() < observations);
        assertTrue(result.randomLinearPredictor("group").length == observations);
        assertTrue(result.smoothTerms().get(0).smoothingParameter() > 0.0);
    }
}
