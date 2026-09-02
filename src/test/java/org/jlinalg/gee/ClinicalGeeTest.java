/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

class ClinicalGeeTest {
    @Test
    void producesAdjustedMeansAndTreatmentContrast() {
        double[] response = {0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0};
        double[][] design = new double[response.length][2];
        int[] cluster = new int[response.length];
        for (int row = 0; row < response.length; row++) {
            design[row][0] = 1.0;
            design[row][1] = row % 2;
            cluster[row] = row / 2;
        }
        GeeResult fit = Gee.fit(response, design, cluster, null,
            GlmFamilies.binomial(), null, null,
            GeeOptions.defaults(), BackendPolicy.CPU);

        MarginalMeanEstimate[] means = ClinicalGee.marginalMeans(fit,
            new double[][] {{1.0, 0.0}, {1.0, 1.0}},
            GlmFamilies.binomial());
        MarginalContrast contrast = ClinicalGee.contrast(fit,
            new double[] {1.0, 1.0}, new double[] {1.0, 0.0});

        assertEquals(2, means.length);
        assertTrue(means[0].mean() > 0.0 && means[0].mean() < 1.0);
        assertTrue(contrast.exponentiatedEstimate() > 0.0);
        assertTrue(contrast.pValue() >= 0.0 && contrast.pValue() <= 1.0);
    }
}
