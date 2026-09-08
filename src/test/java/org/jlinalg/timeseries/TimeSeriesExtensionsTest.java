/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

final class TimeSeriesExtensionsTest {
    @Test
    void integratedAndSparseMissingPathsExposeMetadata() {
        double[] series = new double[80]; for (int index = 1; index < series.length; index++) series[index] = series[index - 1] + Math.sin(index * 0.2) + 0.1;
        DiffuseArima.Result diffuse = DiffuseArima.fit(series, new ArimaOrder(1, 1, 0), ArimaOptions.defaults());
        double[] stationary = new double[40]; for (int index = 0; index < stationary.length; index++) stationary[index] = Math.sin(index * 0.2); stationary[7] = Double.NaN; stationary[23] = Double.NaN;
        SparseMissingSeries.Result sparse = SparseMissingSeries.fit(stationary, new ArimaOrder(1, 0, 0), true, org.jlinalg.compute.BackendPolicy.PREFERRED);
        assertEquals(1, diffuse.diffuseStateCount()); assertEquals(2, sparse.missingCount()); assertTrue(sparse.fit().converged());
    }
}
