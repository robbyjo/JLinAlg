/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.junit.jupiter.api.Test;

/** Checked-in glmmTMB 1.1.14 Laplace reference gates. */
final class ZeroInflatedMixedRReferenceTest {
    @Test
    void groupedZipAndZinbAgreeWithGlmmTmb() throws Exception {
        Data data = data();
        Properties reference = properties();
        RandomEffectTerm count = RandomEffectTerm.randomIntercept(
            "group", data.groups());
        RandomEffectTerm zero = RandomEffectTerm.randomIntercept(
            "group", data.groups());
        ZeroInflatedMixedOptions options = new ZeroInflatedMixedOptions(
            900, 120, 1e-7, 0.35, 1e-8, 1e4,
            1e-5, 1e5, 25.0, null);
        ZeroInflatedMixedResult zip;
        try (SparseZeroInflatedMixedModel.Prepared prepared =
                SparseZeroInflatedMixedModel.preparePoisson(data.zip().length,
                    List.of(count), null, List.of(zero), null, List.of(),
                    options, BackendPolicy.CPU)) {
            zip = prepared.fitWithInference(data.zip(), data.design(), 2,
                data.design(), 2, null, 0, null);
        }
        assertVector(zip.countCoefficients(), reference,
            "zip_count_intercept", "zip_count_x", 1e-4);
        assertVector(zip.zeroCoefficients(), reference,
            "zip_zero_intercept", "zip_zero_x", 1e-4);
        assertEquals(value(reference, "zip_count_sd"),
            Math.sqrt(zip.varianceComponents()[0]), 1e-4);
        assertEquals(value(reference, "zip_zero_sd"),
            Math.sqrt(zip.varianceComponents()[1]), 1e-4);
        assertEquals(value(reference, "zip_log_likelihood"),
            zip.marginalLogLikelihood(), 1e-6);
        double[] standardErrors = zip.standardErrors();
        assertEquals(value(reference, "zip_count_intercept_se"),
            standardErrors[0], 0.025);
        assertEquals(value(reference, "zip_count_x_se"),
            standardErrors[1], 0.02);
        assertEquals(value(reference, "zip_zero_intercept_se"),
            standardErrors[2], 0.05);
        assertEquals(value(reference, "zip_zero_x_se"),
            standardErrors[3], 0.04);

        ZeroInflatedMixedResult zinb =
            SparseZeroInflatedMixedModel.fitNegativeBinomial(data.zinb(),
                data.design(), 2, data.design(), 2, data.intercept(), 1,
                List.of(count), null, List.of(), null, List.of(), null,
                options, BackendPolicy.CPU);
        assertVector(zinb.countCoefficients(), reference,
            "zinb_count_intercept", "zinb_count_x", 1e-4);
        assertVector(zinb.zeroCoefficients(), reference,
            "zinb_zero_intercept", "zinb_zero_x", 1e-4);
        assertEquals(value(reference, "zinb_size"), zinb.sizes()[0], 1e-3);
        assertEquals(value(reference, "zinb_count_sd"),
            Math.sqrt(zinb.varianceComponents()[0]), 1e-4);
        assertEquals(value(reference, "zinb_log_likelihood"),
            zinb.marginalLogLikelihood(), 1e-6);
        System.out.printf("ZI R zip converged=%s LL=%.12f calls=%d; zinb converged=%s LL=%.12f calls=%d%n",
            zip.converged(),zip.marginalLogLikelihood(),zip.objectiveEvaluations(),zinb.converged(),zinb.marginalLogLikelihood(),zinb.objectiveEvaluations());
        assertTrue(zip.converged() && zinb.converged());
    }

    private static void assertVector(
            double[] actual, Properties reference,
            String first, String second, double tolerance) {
        assertEquals(value(reference, first), actual[0], tolerance);
        assertEquals(value(reference, second), actual[1], tolerance);
    }

    private static double value(Properties properties, String name) {
        return Double.parseDouble(properties.getProperty(name));
    }

    private static Properties properties() throws Exception {
        Properties result = new Properties();
        try (InputStream input = ZeroInflatedMixedRReferenceTest.class
                .getResourceAsStream(
                    "/r-reference/zero-inflated-mixed-glmmtmb.properties")) {
            result.load(input);
        }
        return result;
    }

    private static Data data() throws Exception {
        List<Double> zip = new ArrayList<>();
        List<Double> zinb = new ArrayList<>();
        List<Double> x = new ArrayList<>();
        List<String> groups = new ArrayList<>();
        try (InputStream input = ZeroInflatedMixedRReferenceTest.class
                .getResourceAsStream(
                    "/r-reference/zero-inflated-mixed-data.tsv");
             BufferedReader reader = new BufferedReader(new InputStreamReader(
                 input, StandardCharsets.UTF_8))) {
            reader.readLine();
            for (String line = reader.readLine(); line != null;
                    line = reader.readLine()) {
                String[] fields = line.split("\t", -1);
                zip.add(Double.valueOf(fields[0]));
                zinb.add(Double.valueOf(fields[1]));
                x.add(Double.valueOf(fields[2]));
                groups.add(fields[3]);
            }
        }
        int rows = zip.size();
        double[] design = new double[rows * 2];
        double[] intercept = new double[rows];
        for (int row = 0; row < rows; row++) {
            design[2 * row] = 1.0;
            design[2 * row + 1] = x.get(row);
            intercept[row] = 1.0;
        }
        return new Data(toArray(zip), toArray(zinb), design, intercept,
            List.copyOf(groups));
    }

    private static double[] toArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int index = 0; index < result.length; index++)
            result[index] = values.get(index);
        return result;
    }

    private record Data(
            double[] zip, double[] zinb, double[] design,
            double[] intercept, List<String> groups) { }
}
