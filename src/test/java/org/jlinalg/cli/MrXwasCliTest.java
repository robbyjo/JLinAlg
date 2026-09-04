/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MrXwasCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void scansLongFormatExposureAndMultiplePhenotypeFiles() throws Exception {
        Path exposure = temporaryDirectory.resolve("clumped.tsv");
        Path outcome = temporaryDirectory.resolve("outcomes.tsv");
        Path result = temporaryDirectory.resolve("mr-results.tsv");
        Path fdr = temporaryDirectory.resolve("mr-all-pairs.tsv");
        Files.writeString(exposure, exposureTable(), StandardCharsets.UTF_8);
        Files.writeString(outcome, outcomeTable(), StandardCharsets.UTF_8);
        ByteArrayOutputStream standard = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();

        int status = JLinAlgCli.run(new String[] {"mr-xwas",
            "--exposure", exposure.toString(), "--outcome", outcome.toString(),
            "--output", result.toString(), "--p-threshold", "1e-4",
            "--fdr-output", fdr.toString(),
            "--threads", "4", "--pair-block-size", "3",
            "--bootstrap-replicates", "50", "--seed", "42"
        }, new PrintStream(standard), new PrintStream(error));

        assertEquals(0, status, error.toString(StandardCharsets.UTF_8));
        List<String> lines = Files.readAllLines(result);
        assertEquals(3, lines.size());
        assertTrue(lines.get(0).contains("negative_log10_p_value"));
        assertTrue(lines.get(1).startsWith(
            "GENE1\tGENE1\tcardiovascular\tCAD\tCoronary artery disease\t"));
        assertTrue(lines.get(2).startsWith(
            "PROT1\tPROT1\tkidney\tCKD\tChronic kidney disease\t"));
        Path failures = temporaryDirectory.resolve("mr-results.failures.tsv");
        assertEquals(1, Files.readAllLines(failures).size());
        List<String> adjusted = Files.readAllLines(fdr);
        assertEquals(5, adjusted.size());
        assertTrue(adjusted.get(0).endsWith("\tfdr_bh"));
        assertTrue(adjusted.get(1).startsWith(
            "GENE1\tGENE1\tcardiovascular\tCAD\t"));
        assertTrue(adjusted.get(2).startsWith(
            "GENE1\tGENE1\tkidney\tCKD\t"));
        assertTrue(adjusted.get(3).startsWith(
            "PROT1\tPROT1\tcardiovascular\tCAD\t"));
        assertTrue(adjusted.get(4).startsWith(
            "PROT1\tPROT1\tkidney\tCKD\t"));
        double[] pValues = adjusted.stream().skip(1)
            .mapToDouble(line -> Double.parseDouble(line.split("\t")[10]))
            .toArray();
        double[] qValues = adjusted.stream().skip(1)
            .mapToDouble(line -> Double.parseDouble(line.split("\t")[20]))
            .toArray();
        double[] expected = bh(pValues);
        for (int index = 0; index < expected.length; index++)
            assertEquals(expected[index], qValues[index], 1e-14);
        String report = standard.toString(StandardCharsets.UTF_8);
        assertTrue(report.contains("6 exposure-outcome pairs (2 x 3)"));
        assertTrue(report.contains("Retained 2; below threshold 2"));
        assertTrue(report.contains("fewer than 3 harmonized instruments 2"));
        assertTrue(report.contains("BH-adjusted 4 successfully screened pairs"));
    }

    @Test
    void rejectsAmbiguousThresholdScales() throws Exception {
        Path exposure = temporaryDirectory.resolve("clumped.tsv");
        Path outcome = temporaryDirectory.resolve("outcomes.tsv");
        Files.writeString(exposure, exposureTable(), StandardCharsets.UTF_8);
        Files.writeString(outcome, outcomeTable(), StandardCharsets.UTF_8);
        ByteArrayOutputStream messages = new ByteArrayOutputStream();

        int status = JLinAlgCli.run(new String[] {"mr-xwas",
            "--exposure", exposure.toString(), "--outcome", outcome.toString(),
            "--output", temporaryDirectory.resolve("out.tsv").toString(),
            "--p-threshold", "0.05", "--negative-log10-p-threshold", "2"
        }, new PrintStream(messages), new PrintStream(messages));

        assertEquals(2, status);
        assertTrue(messages.toString(StandardCharsets.UTF_8)
            .contains("specify exactly one p-value threshold scale"));
    }

    private static String exposureTable() {
        StringBuilder result = new StringBuilder(
            "Phenotype\tSNP\tbeta\tse\teaf\teffect_allele\tother_allele\tpval\tgene\n");
        appendExposure(result, "GENE1", "a", new double[] {0.10, 0.15, 0.20});
        appendExposure(result, "PROT1", "b", new double[] {0.08, 0.12, 0.17});
        return result.toString();
    }

    private static void appendExposure(StringBuilder result, String id,
            String prefix, double[] effects) {
        for (int index = 0; index < effects.length; index++)
            result.append(id).append('\t').append(prefix).append(index + 1)
                .append('\t').append(effects[index])
                .append("\t0.01\t0.2\tA\tC\t1e-9\t\n");
    }

    private static String outcomeTable() {
        StringBuilder result = new StringBuilder(
            "id.outcome\toutcome\tcategory\tSNP\tbeta\tse\teaf\teffect_allele\tother_allele\tpval\n");
        appendOutcome(result, "CAD", "Coronary artery disease",
            "cardiovascular", 0.5, 0.0);
        appendOutcome(result, "CKD", "Chronic kidney disease",
            "kidney", 0.0, -0.4);
        result.append("COPD\tChronic obstructive pulmonary disease\tlung\ta1")
            .append("\t0.01\t0.1\t0.2\tA\tC\t0.9\n")
            .append("COPD\tChronic obstructive pulmonary disease\tlung\tb1")
            .append("\t0.01\t0.1\t0.2\tA\tC\t0.9\n");
        return result.toString();
    }

    private static void appendOutcome(StringBuilder result, String id,
            String label, String category, double firstEffect,
            double secondEffect) {
        double[] first = {0.10, 0.15, 0.20};
        double[] second = {0.08, 0.12, 0.17};
        for (int index = 0; index < first.length; index++)
            outcomeRow(result, id, label, category, "a" + (index + 1),
                firstEffect * first[index], firstEffect == 0.0 ? 0.1 : 0.01);
        for (int index = 0; index < second.length; index++)
            outcomeRow(result, id, label, category, "b" + (index + 1),
                secondEffect * second[index], secondEffect == 0.0 ? 0.1 : 0.01);
    }

    private static void outcomeRow(StringBuilder result, String id,
            String label, String category, String variant, double effect,
            double standardError) {
        result.append(id).append('\t').append(label).append('\t')
            .append(category).append('\t').append(variant).append('\t')
            .append(effect).append('\t').append(standardError)
            .append("\t0.2\tA\tC\t0.5\n");
    }

    private static double[] bh(double[] pValues) {
        Integer[] order = new Integer[pValues.length];
        for (int index = 0; index < order.length; index++) order[index] = index;
        java.util.Arrays.sort(order,
            java.util.Comparator.comparingDouble(index -> pValues[index]));
        double[] result = new double[pValues.length];
        double running = 1.0;
        for (int rank = order.length; rank >= 1; rank--) {
            int index = order[rank - 1];
            running = Math.min(running,
                pValues[index] * order.length / rank);
            result[index] = Math.min(1.0, running);
        }
        return result;
    }
}
