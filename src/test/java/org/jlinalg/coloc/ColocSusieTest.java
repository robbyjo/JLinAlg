/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.susie.Susie;
import org.jlinalg.susie.SusieOptions;
import org.jlinalg.susie.SusieResult;
import org.junit.jupiter.api.Test;

class ColocSusieTest {
    @Test
    void matchesColocSusieForBundledD1D2Example() throws IOException {
        Fixture fixture = fixture();
        ColocSusieResult result = ColocSusie.analyze(
            fixture.input("D1.L1"), fixture.input("D2.L1"));

        assertEquals(500, result.commonVariants().size());
        assertEquals(1, result.signalPairs().size());
        assertEquals(0, result.skippedSignalPairs());
        ColocSignalPair pair = result.signalPairs().get(0);
        assertEquals("s105", pair.trait1LeadVariant());
        assertEquals("s105", pair.trait2LeadVariant());
        assertArrayEquals(new double[] {
            2.139563012548298e-19,
            1.743392886839919e-10,
            4.334120750809625e-12,
            0.0015346667878764,
            0.9984653330334475
        }, pair.hypothesisPosteriors(), 2e-12);

        double[] shared = result.sharedVariantPosterior(0);
        assertEquals(0.9897938695610928, shared[104], 2e-12);
        assertEquals(0.01015156621174276, shared[102], 2e-12);
        assertEquals(1.0, sum(shared), 2e-15);
    }

    @Test
    void matchesColocSusieForBundledMultiSignalD3D4Example()
            throws IOException {
        Fixture fixture = fixture();
        ColocSusieInput d3 = fixture.input("D3.L1", "D3.L2");
        ColocSusieResult result = ColocSusie.analyze(
            d3, fixture.input("D4.L1"));

        assertEquals(2, result.signalPairs().size());
        ColocSignalPair first = result.signalPairs().get(0);
        ColocSignalPair second = result.signalPairs().get(1);
        assertEquals("s105", first.trait1LeadVariant());
        assertEquals("s89", second.trait1LeadVariant());
        assertEquals(0, first.trait1EffectIndex());
        assertEquals(1, second.trait1EffectIndex());
        assertArrayEquals(new double[] {
            3.079007921550071e-14,
            6.507290798949911e-7,
            1.342029586955264e-10,
            0.0008379729004194668,
            0.99916137623626555
        }, first.hypothesisPosteriors(), 2e-12);
        assertArrayEquals(new double[] {
            1.422896088646072e-6,
            2.209787034466888e-4,
            0.006201895866395195,
            0.9631063075065081,
            0.03046939502756055
        }, second.hypothesisPosteriors(), 2e-12);
        assertEquals(0.9983776660042696,
            result.sharedVariantPosterior(0)[104], 2e-12);
        assertEquals(0.6558882769068253,
            result.sharedVariantPosterior(1)[104], 2e-12);
    }

    @Test
    void alignsVariantsAndTrimsPairsWithoutPosteriorCoverage() {
        ColocSusieInput first = new ColocSusieInput(
            List.of("a", "b", "c"), new double[][] {{10.0, 0.0, 0.0}});
        ColocSusieInput second = new ColocSusieInput(
            List.of("b", "c", "d"), new double[][] {{0.0, 0.0, 10.0}});

        ColocSusieResult result = ColocSusie.analyze(first, second);

        assertEquals(List.of("b", "c"), result.commonVariants());
        assertTrue(result.signalPairs().isEmpty());
        assertEquals(1, result.skippedSignalPairs());
    }

    @Test
    void matchesColocWeightedPriorEquations() {
        List<String> variants = List.of("a", "b", "c");
        ColocSusieInput first = new ColocSusieInput(
            variants, new double[][] {{2.0, -1.0, 0.5}});
        ColocSusieInput second = new ColocSusieInput(
            variants, new double[][] {{1.0, 0.2, -0.5}});
        ColocOptions options = new ColocOptions(
            1e-4, 1e-4, 5e-6, 0.5, true,
            new double[] {1.0, 2.0, 3.0},
            new double[] {3.0, 1.0, 2.0});

        ColocSusieResult result = ColocSusie.analyze(first, second, options);

        assertArrayEquals(new double[] {
            0.99873446695617885,
            0.0006527218518897913,
            0.0005287954220899698,
            1.7791798369307635e-7,
            8.38378518577422e-5
        }, result.signalPairs().get(0).hypothesisPosteriors(), 2e-15);
        assertArrayEquals(new double[] {
            0.9326984912884172,
            0.020865185159203038,
            0.04643632355237998
        }, result.sharedVariantPosterior(0), 2e-15);
    }

    @Test
    void impossibleSharedConfigurationFallsBackToNull() {
        List<String> variants = List.of("a", "b");
        ColocSusieInput first = new ColocSusieInput(variants,
            new double[][] {{0.0, Double.NEGATIVE_INFINITY}});
        ColocSusieInput second = new ColocSusieInput(variants,
            new double[][] {{Double.NEGATIVE_INFINITY, 0.0}});
        ColocOptions options = new ColocOptions(
            1e-4, 1e-4, 5e-6, 0.5, false, null, null);

        ColocSusieResult result = ColocSusie.analyze(first, second, options);

        assertArrayEquals(new double[] {1.0, 0.0, 0.0, 0.0, 0.0},
            result.signalPairs().get(0).hypothesisPosteriors());
        assertArrayEquals(new double[] {0.0, 0.0},
            result.sharedVariantPosterior(0));
    }

    @Test
    void consumesJlinAlgSusieCredibleEffectsDirectly() {
        double[] z = {12.0, 0.2, -0.1, 0.3};
        double[][] ld = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0, 1.0}
        };
        SusieResult fit = Susie.fitSummary(z, ld, 1_000,
            List.of("a", "b", "c", "d"),
            new SusieOptions(1, 100, 1e-8, 0.2, false, 0.95, 0.0),
            BackendPolicy.CPU);

        ColocSusieResult result = ColocSusie.analyze(fit, fit);

        assertEquals(1, result.signalPairs().size());
        assertEquals("a", result.signalPairs().get(0).trait1LeadVariant());
        assertTrue(result.signalPairs().get(0).posteriorH4() > 0.99);
        assertTrue(result.sharedVariantPosterior(0)[0] > 0.999999);
    }

    private static double sum(double[] values) {
        double result = 0.0;
        for (double value : values) result += value;
        return result;
    }

    private static Fixture fixture() throws IOException {
        String path = "/r-reference/coloc-susie-example.tsv";
        InputStream stream = ColocSusieTest.class.getResourceAsStream(path);
        if (stream == null) throw new IOException("missing resource " + path);
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(
                stream, StandardCharsets.UTF_8))) {
            String header;
            do {
                header = reader.readLine();
            } while (header != null && header.startsWith("#"));
            if (header == null) throw new IOException("missing fixture header");
            String[] columns = header.split("\\t");
            List<String> variants = new ArrayList<>();
            Map<String, List<Double>> values = new LinkedHashMap<>();
            for (int column = 1; column < columns.length; column++) {
                values.put(columns[column], new ArrayList<>());
            }
            for (String line = reader.readLine(); line != null;
                    line = reader.readLine()) {
                String[] fields = line.split("\\t");
                variants.add(fields[0]);
                for (int column = 1; column < fields.length; column++) {
                    values.get(columns[column]).add(
                        Double.parseDouble(fields[column]));
                }
            }
            Map<String, double[]> arrays = new LinkedHashMap<>();
            values.forEach((name, source) -> {
                double[] target = new double[source.size()];
                for (int index = 0; index < target.length; index++) {
                    target[index] = source.get(index);
                }
                arrays.put(name, target);
            });
            return new Fixture(List.copyOf(variants), arrays);
        }
    }

    private record Fixture(List<String> variants, Map<String, double[]> values) {
        ColocSusieInput input(String... signals) {
            double[][] matrix = new double[signals.length][];
            for (int index = 0; index < signals.length; index++) {
                matrix[index] = values.get(signals[index]);
            }
            return new ColocSusieInput(variants, matrix);
        }
    }
}
