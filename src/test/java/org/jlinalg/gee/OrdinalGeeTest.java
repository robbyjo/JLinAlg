/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class OrdinalGeeTest {
    @Test
    void proportionalOddsFitReturnsThresholdsAndSlope() {
        int clusters = 30;
        int visits = 2;
        int[] response = new int[clusters * visits];
        double[][] covariates = new double[response.length][1];
        int[] id = new int[response.length];
        int[] repeated = new int[response.length];
        for (int cluster = 0; cluster < clusters; cluster++) {
            for (int visit = 0; visit < visits; visit++) {
                int row = cluster * visits + visit;
                response[row] = (cluster + visit) % 3;
                covariates[row][0] = (cluster % 4) - 1.5;
                id[row] = cluster;
                repeated[row] = visit;
            }
        }

        OrdinalGeeResult result = OrdinalGee.fit(response, covariates,
            id, repeated, 3, GeeOptions.defaults(), BackendPolicy.CPU);

        assertEquals(2, result.thresholds().length);
        assertEquals(1, result.coefficients().length);
        assertTrue(Arrays.stream(result.fit().covariance()).allMatch(Double::isFinite));
    }
}
