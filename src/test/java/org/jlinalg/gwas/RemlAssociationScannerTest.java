/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gwas;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.Ols;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class RemlAssociationScannerTest {
    @Test
    void batchedP3dEffectsMatchPerMarkerOlsForIdentityCovariance() {
        double[] response = {1, 2, 2, 4, 5, 7, 7, 9};
        double[][] covariates = {
            {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
        };
        double[][] markers = {
            {0, 0}, {0, 1}, {1, 0}, {1, 1},
            {2, 0}, {2, 1}, {3, 0}, {3, 1}
        };
        RemlAssociationScanner scanner = RemlAssociationScanner.prepare(
            response, covariates,
            List.of(VarianceComponent.identity("residual", response.length)),
            RemlOptions.builder().initialVariances(2.0).build(),
            BackendPolicy.CPU);

        AssociationScanResult scan = scanner.scan(markers, List.of("g1", "g2"),
            new AssociationScanOptions(
                1, GenotypeMissingPolicy.MEAN_IMPUTE, 2));

        for (int marker = 0; marker < 2; marker++) {
            double[][] full = new double[response.length][2];
            for (int row = 0; row < response.length; row++) {
                full[row][0] = 1.0;
                full[row][1] = markers[row][marker];
            }
            assertEquals(Ols.fit(response, full).coefficients()[1],
                scan.beta()[marker], 1e-10);
            assertTrue(Double.isFinite(scan.pValues()[marker]));
        }
    }
}
