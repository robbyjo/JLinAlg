/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.io.InputStream;
import java.nio.file.Path;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mr.MultivariateMendelianRandomization;
import org.junit.jupiter.api.Test;

final class MultivariateMrReferenceTest {
    @Test void jointIvwMatchesIndependentBaseRMatrixCalculation()
            throws Exception {
        Path root = Path.of("src/test/resources/r-reference");
        var instruments = MrWideTable.read(
            root.resolve("multivariate-mr-input.tsv"),
            List.of("x1", "x2"), List.of("y1", "y2"));
        var result = MultivariateMendelianRandomization.fit(instruments,
            List.of("x1", "x2"), List.of("y1", "y2"),
            MrWideTable.readMatrix(
                root.resolve("multivariate-mr-correlation.tsv")),
            false, BackendPolicy.CPU);
        Properties reference = new Properties();
        try (InputStream stream = getClass().getResourceAsStream(
                "/r-reference/multivariate-mr-reference.properties")) {
            assertNotNull(stream);
            reference.load(stream);
        }
        assertValues(reference, "beta", result.beta(), 2e-13);
        assertValues(reference, "se", result.standardErrors(), 2e-13);
        assertValues(reference, "covariance", result.covariance(), 3e-13);
        assertEquals(number(reference, "q"), result.cochranQ(), 3e-13);
        assertEquals((int) number(reference, "q_df"),
            result.heterogeneityDegreesOfFreedom());
        assertEquals(number(reference, "q_p"),
            result.heterogeneityPValue(), 2e-13);
        assertEquals(number(reference, "overall"),
            result.overallTest().chiSquare(), 3e-13);
        assertEquals(number(reference, "overall_p"),
            result.overallTest().pValue(), 2e-13);
        assertEquals(number(reference, "exposure_x1"),
            result.exposureJointTests().get(0).chiSquare(), 3e-13);
        assertEquals(number(reference, "exposure_x2"),
            result.exposureJointTests().get(1).chiSquare(), 3e-13);
        assertEquals(number(reference, "outcome_y1"),
            result.outcomeJointTests().get(0).chiSquare(), 3e-13);
        assertEquals(number(reference, "outcome_y2"),
            result.outcomeJointTests().get(1).chiSquare(), 3e-13);
    }

    private static void assertValues(Properties reference, String prefix,
            double[] actual, double tolerance) {
        for (int index = 0; index < actual.length; index++)
            assertEquals(number(reference, prefix + "." + index), actual[index],
                tolerance, prefix + "[" + index + "]");
    }

    private static double number(Properties reference, String name) {
        return Double.parseDouble(reference.getProperty(name));
    }
}
