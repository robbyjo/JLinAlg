/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

final class MgcvReferenceTest {
    @Test
    void gaussianPsplineTracksCommittedMgcvFit() throws IOException {
        Properties reference = new Properties();
        try (InputStream input = getClass().getResourceAsStream(
                "/r-reference/mgcv-gaussian-reml.properties")) {
            if (input == null) throw new IOException("mgcv fixture is missing");
            reference.load(input);
        }
        int observations = 80;
        double[] x = new double[observations];
        double[] y = new double[observations];
        double[][] intercept = new double[observations][1];
        for (int row = 0; row < observations; row++) {
            x[row] = row / 79.0;
            y[row] = 1.5 + Math.sin(2.0 * Math.PI * x[row])
                + 0.08 * Math.sin(17.0 * (row + 1));
            intercept[row][0] = 1.0;
        }
        GamResult fit = Gam.fitGaussian(y, intercept,
            List.of(PSplineTerm.of("s(x)", x, 10)),
            RemlOptions.defaults(), BackendPolicy.CPU);
        double[] mgcv = parse(reference.getProperty("fitted"));
        assertEquals(observations, mgcv.length);
        double squaredDifference = 0.0;
        for (int row = 0; row < observations; row++) {
            double difference = fit.fittedValues()[row] - mgcv[row];
            squaredDifference += difference * difference;
        }
        double rmse = Math.sqrt(squaredDifference / observations);
        assertTrue(rmse < 0.035, "JLinAlg/mgcv fitted RMSE=" + rmse);
        // mgcv's s.table EDF excludes the model intercept; JLinAlg totals it.
        assertEquals(1.0 + Double.parseDouble(reference.getProperty("edf")),
            fit.totalEffectiveDegreesOfFreedom(), 1.0);
    }

    private static double[] parse(String values) {
        String[] fields = values.split(",");
        double[] result = new double[fields.length];
        for (int index = 0; index < fields.length; index++) {
            result[index] = Double.parseDouble(fields[index]);
        }
        return result;
    }
}
