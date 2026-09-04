/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class SusieRReferenceTest {
    private static final byte[] MAGIC = "JLSUSIE1".getBytes(StandardCharsets.US_ASCII);

    @Test
    void matchesSusieROnN3FineMappingVignetteData() throws IOException {
        Data data = readData();
        Reference reference = readReference();
        List<String> names = new ArrayList<>(data.variables());
        for (int index = 0; index < data.variables(); index++) {
            names.add("variable" + (index + 1));
        }

        SusieResult result = Susie.fit(data.response(), data.design(), names,
            new SusieOptions(10, 200, 1e-6, 0.2, true, 0.95, 0.5),
            BackendPolicy.CPU);

        assertTrue(result.converged());
        assertEquals(reference.iterations(), result.iterations());
        assertEquals(reference.residualVariance(), result.residualVariance(), 2e-10);
        assertEquals(reference.intercept(), result.intercept(), 2e-10);
        assertEquals(reference.objective(), result.objective(), 2e-8);
        assertEquals(3, result.credibleSets().size());
        assertMaximumError(reference.pip(), result.pip(), 2e-10);
        assertMaximumError(reference.coefficient(), result.posteriorMean(), 2e-10);
    }

    private static void assertMaximumError(
            double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        double maximum = 0.0;
        for (int index = 0; index < expected.length; index++) {
            maximum = Math.max(maximum, Math.abs(expected[index] - actual[index]));
        }
        assertTrue(maximum <= tolerance,
            "maximum absolute error " + maximum + " exceeds " + tolerance);
    }

    private static Data readData() throws IOException {
        try (DataInputStream input = new DataInputStream(new BufferedInputStream(
                new GZIPInputStream(resource("/susie/N3finemapping.bin.gz"))))) {
            assertTrue(Arrays.equals(MAGIC, input.readNBytes(MAGIC.length)));
            int observations = input.readInt();
            int variables = input.readInt();
            double[][] design = new double[observations][variables];
            for (double[] row : design) {
                for (int column = 0; column < variables; column++) {
                    row[column] = input.readDouble();
                }
            }
            double[] response = new double[observations];
            for (int row = 0; row < observations; row++) {
                response[row] = input.readDouble();
            }
            for (int column = 0; column < variables; column++) input.readDouble();
            assertEquals(-1, input.read());
            return new Data(design, response);
        }
    }

    private static Reference readReference() throws IOException {
        double residualVariance = Double.NaN;
        double intercept = Double.NaN;
        double objective = Double.NaN;
        int iterations = -1;
        List<Double> pip = new ArrayList<>();
        List<Double> coefficient = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(
                resource("/susie/N3finemapping-reference.tsv"),
                StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("# residual_variance=")) {
                    residualVariance = value(line);
                } else if (line.startsWith("# intercept=")) {
                    intercept = value(line);
                } else if (line.startsWith("# iterations=")) {
                    iterations = (int) value(line);
                } else if (line.startsWith("# objective=")) {
                    objective = value(line);
                } else if (!line.startsWith("#") && !line.startsWith("index")) {
                    String[] fields = line.split("\\t", -1);
                    assertEquals(pip.size() + 1, Integer.parseInt(fields[0]));
                    pip.add(Double.parseDouble(fields[1]));
                    coefficient.add(Double.parseDouble(fields[2]));
                }
            }
        }
        return new Reference(toArray(pip), toArray(coefficient),
            residualVariance, intercept, iterations, objective);
    }

    private static double value(String line) {
        return Double.parseDouble(line.substring(line.indexOf('=') + 1));
    }

    private static double[] toArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int index = 0; index < result.length; index++) result[index] = values.get(index);
        return result;
    }

    private static InputStream resource(String path) {
        InputStream result = SusieRReferenceTest.class.getResourceAsStream(path);
        if (result == null) throw new IllegalStateException("missing resource " + path);
        return result;
    }

    private record Data(double[][] design, double[] response) {
        int variables() { return design[0].length; }
    }

    private record Reference(
            double[] pip, double[] coefficient, double residualVariance,
            double intercept, int iterations, double objective) { }
}
