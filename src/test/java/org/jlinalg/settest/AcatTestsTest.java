/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;

class AcatTestsTest {
    @Test void genericCombinationMatchesReferenceAndStabilizesTails() {
        Acat.Result result = Acat.combine(
            new double[] {.02, .0004, .2, .1, .8});
        assertEquals(0.0019534044057701767, result.pValue(), 2e-18);
        assertEquals(5, result.components());

        Acat.Result tail = Acat.combine(
            new double[] {1e-300, .2}, new double[] {3, 1});
        assertEquals(1e-300 / .75, tail.pValue(), 1e-314);
        assertThrows(IllegalArgumentException.class,
            () -> Acat.combine(new double[] {0, 1}));
        assertThrows(IllegalArgumentException.class,
            () -> Acat.combine(new double[] {.1, .2}, new double[] {0, 0}));
    }

    @Test void acatVAndCanonicalAcatOMatchIndependentRFixtures()
            throws Exception {
        Path fixture = Path.of("src/test/resources/acat/reference.tsv");
        var lines = Files.readAllLines(fixture);
        String[] header = lines.get(0).split("\t");
        Map<String, Integer> columns = java.util.stream.IntStream
            .range(0, header.length).boxed()
            .collect(Collectors.toMap(index -> header[index], index -> index));
        for (String line : lines.subList(1, lines.size())) {
            String[] fields = line.split("\t");
            String id = fields[columns.get("id")];
            double[] scores = numbers(fields[columns.get("u")]);
            double[] information = numbers(fields[columns.get("v")]);
            double[] mafs = numbers(fields[columns.get("maf")]);
            double[] macs = numbers(fields[columns.get("mac")]);
            SetTestScoreState state = new SetTestScoreState(
                scores, information, scores.length);

            AcatVResult rare = SummarySetTests.acatV(
                id, state, mafs, macs, 1, 25, 10);
            AcatVResult equal = SummarySetTests.acatV(
                id, state, mafs, macs, 1, 1, 10);
            assertClose(fields, columns, "acat_v_1_25", rare.pValue(), id);
            assertClose(fields, columns, "acat_v_1_1", equal.pValue(), id);
            assertEquals((int) Arrays.stream(macs).filter(value -> value <= 10).count(),
                rare.collapsedVariants(), id);
            assertEquals(rare.collapsedVariants() > 0 ? 1 : 0,
                rare.components().stream()
                .filter(component -> component.name().equals("ultra-rare-burden"))
                .count(), id);

            AcatOResult omnibus = SummarySetTests.acatO(
                id, state, mafs, macs, 10);
            assertClose(fields, columns, "acat_o", omnibus.pValue(), id, 5e-5, 1e-5);
            assertEquals(6, omnibus.components().size(), id);
            Map<String, Double> componentP = omnibus.components().stream()
                .collect(Collectors.toMap(
                    AcatOResult.Component::name, AcatOResult.Component::pValue));
            assertClose(fields, columns, "skat_1_25",
                componentP.get("skat-beta-1-25"), id, 5e-5, 7e-5);
            assertClose(fields, columns, "skat_1_1",
                componentP.get("skat-beta-1-1"), id, 5e-5, 7e-5);
            assertClose(fields, columns, "burden_1_25",
                componentP.get("burden-beta-1-25"), id);
            assertClose(fields, columns, "burden_1_1",
                componentP.get("burden-beta-1-1"), id);
        }
    }

    @Test void invalidMomentsAndDimensionsFailExplicitly() {
        SetTestScoreState state = new SetTestScoreState(
            new double[] {1, 2}, new double[] {1, 2, 2, 1}, 2);
        assertThrows(IllegalArgumentException.class,
            () -> AcatTests.acatV("g", state,
                new double[] {.01, .02}, new double[] {5, 20}));
        assertThrows(IllegalArgumentException.class,
            () -> AcatTests.acatV("g", state,
                new double[] {.01}, new double[] {5}));
        assertThrows(IllegalArgumentException.class,
            () -> AcatTests.acatO("g", state,
                new double[] {.01, .02}, new double[] {5, 20}, -1));
    }

    @Test void mixedAcatVMatchesReferenceExactOneFiltering() {
        SetTestScoreState state = new SetTestScoreState(
            new double[] {0, 1}, new double[] {1, 0, 0, 1}, 2);

        AcatVResult result = AcatTests.acatV("g", state,
            new double[] {.001, .02}, new double[] {2, 40}, 1, 1, 10);

        assertEquals(1, result.components().size());
        assertEquals("variant-2", result.components().get(0).name());
        assertEquals(2 * jdistlib.Normal.cumulative(-1, 0, 1, true, false),
            result.pValue(), 1e-15);
    }

    private static double[] numbers(String value) {
        return Arrays.stream(value.split(","))
            .map(String::trim).mapToDouble(Double::parseDouble).toArray();
    }

    private static void assertClose(
            String[] fields, Map<String, Integer> columns, String column,
            double actual, String id) {
        assertClose(fields, columns, column, actual, id, 2e-10);
    }

    private static void assertClose(
            String[] fields, Map<String, Integer> columns, String column,
            double actual, String id, double relativeTolerance) {
        assertClose(fields, columns, column, actual, id, relativeTolerance, 2e-12);
    }

    private static void assertClose(
            String[] fields, Map<String, Integer> columns, String column,
            double actual, String id, double relativeTolerance,
            double absoluteTolerance) {
        double expected = Double.parseDouble(fields[columns.get(column)]);
        assertEquals(expected, actual,
            Math.max(relativeTolerance * expected, absoluteTolerance),
            id + " " + column);
    }
}
