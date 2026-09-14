/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Properties;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class InstrumentalVariableRReferenceTest {
    private static final String ROOT = "/r-reference/instrumental-variables-";
    private static final Map<String, double[]> REFERENCE = reference();

    @Test
    void coefficientsCovariancesAndStrengthMatchIndependentBaseR() {
        Data data = data();
        for (InstrumentalVariableCovariance estimator
                : InstrumentalVariableCovariance.values()) {
            InstrumentalVariableResult fit = InstrumentalVariableRegression.fit(
                data.response(), data.exogenous(), data.endogenous(),
                data.instruments(), estimator.requiresClusters()
                    ? data.clusters() : null,
                new InstrumentalVariableOptions(estimator, 0.95),
                BackendPolicy.CPU);
            String prefix = estimator.name().toLowerCase();
            assertArrayEquals(REFERENCE.get("beta"), fit.coefficients(),
                2e-11, estimator.name());
            assertArrayEquals(REFERENCE.get(prefix + "_cov"),
                fit.covariance(), 3e-11, estimator.name());
            assertEquals(REFERENCE.get("rss")[0],
                fit.residualSumOfSquares(), 2e-11);
            assertEquals(REFERENCE.get("residual_variance")[0],
                fit.residualVariance(), 2e-13);
            for (int column = 0; column < 2; column++) {
                double[] expected = REFERENCE.get(prefix
                    + "_first_stage_" + column);
                InstrumentStrengthDiagnostic actual = fit.strengthDiagnostics()
                    .get(column);
                assertClose(expected[0], actual.partialRSquared(), 3e-11);
                assertClose(expected[1], actual.classicalFStatistic(), 3e-9);
                assertClose(expected[2], actual.waldStatistic(), 2e-8);
                assertClose(expected[3], actual.effectiveFStatistic(), 1e-8);
                assertClose(expected[4], actual.pValue(), 2e-10);
            }
        }
    }

    @Test
    void projectedValuesAndConditionDiagnosticMatchRQrAndSvd() {
        Data data = data();
        InstrumentalVariableResult fit = InstrumentalVariableRegression.fit(
            data.response(), data.exogenous(), data.endogenous(),
            data.instruments(), InstrumentalVariableOptions.defaults(),
            BackendPolicy.CPU);
        int[] selected = {0, 1, 17, 78, 179};
        double[][] fitted = fit.instrumentedEndogenous();
        double[] actual = new double[selected.length * 2];
        for (int index = 0; index < selected.length; index++) {
            actual[index * 2] = fitted[selected[index]][0];
            actual[index * 2 + 1] = fitted[selected[index]][1];
        }
        assertArrayEquals(REFERENCE.get("instrumented_selected"), actual,
            2e-12);
        assertEquals(REFERENCE.get("projected_condition")[0],
            fit.projectedDesignConditionNumber(), 2e-12);
    }

    @Test
    void fixtureRecordsIndependentRRuntimeAndImplementation() throws IOException {
        Properties properties = new Properties();
        try (var stream = InstrumentalVariableRReferenceTest.class
                .getResourceAsStream(ROOT + "reference.properties")) {
            properties.load(Objects.requireNonNull(stream));
        }
        assertTrue(properties.getProperty("runtime").startsWith("R version"));
        assertTrue(properties.getProperty("implementation")
            .contains("explicit IV/GMM sandwiches"));
    }

    private static void assertClose(
            double expected, double actual, double tolerance) {
        assertEquals(expected, actual,
            tolerance * Math.max(1.0, Math.abs(expected)));
    }

    private static Data data() {
        List<double[]> rows = new ArrayList<>();
        try (var stream = InstrumentalVariableRReferenceTest.class
                .getResourceAsStream(ROOT + "data.csv");
             var reader = new BufferedReader(new InputStreamReader(
                 Objects.requireNonNull(stream), StandardCharsets.UTF_8))) {
            reader.readLine();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split(",");
                double[] values = new double[fields.length];
                for (int index = 0; index < fields.length; index++) {
                    values[index] = Double.parseDouble(fields[index]);
                }
                rows.add(values);
            }
        } catch (IOException exception) {
            throw new UncheckedIOException(exception);
        }
        int observations = rows.size();
        double[] response = new double[observations];
        double[][] exogenous = new double[observations][2];
        double[][] endogenous = new double[observations][2];
        double[][] instruments = new double[observations][3];
        int[] clusters = new int[observations];
        for (int row = 0; row < observations; row++) {
            double[] values = rows.get(row);
            clusters[row] = (int) values[0];
            exogenous[row][0] = 1.0;
            exogenous[row][1] = values[1];
            instruments[row][0] = values[2];
            instruments[row][1] = values[3];
            instruments[row][2] = values[4];
            endogenous[row][0] = values[5];
            endogenous[row][1] = values[6];
            response[row] = values[7];
        }
        return new Data(response, exogenous, endogenous, instruments, clusters);
    }

    private static Map<String, double[]> reference() {
        Map<String, List<Double>> values = new HashMap<>();
        try (var stream = InstrumentalVariableRReferenceTest.class
                .getResourceAsStream(ROOT + "reference.csv");
             var reader = new BufferedReader(new InputStreamReader(
                 Objects.requireNonNull(stream), StandardCharsets.UTF_8))) {
            reader.readLine();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split(",");
                values.computeIfAbsent(fields[0], ignored -> new ArrayList<>())
                    .add(Double.valueOf(fields[2]));
            }
        } catch (IOException exception) {
            throw new UncheckedIOException(exception);
        }
        Map<String, double[]> result = new HashMap<>();
        values.forEach((key, entries) -> result.put(key,
            entries.stream().mapToDouble(Double::doubleValue).toArray()));
        return result;
    }

    private record Data(
            double[] response,
            double[][] exogenous,
            double[][] endogenous,
            double[][] instruments,
            int[] clusters) { }
}
