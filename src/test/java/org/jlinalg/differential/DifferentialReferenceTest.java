/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.differential;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

final class DifferentialReferenceTest {
    private static final double[][] DESIGN = {
        {1, 0}, {1, 0}, {1, 0}, {1, 0},
        {1, 1}, {1, 1}, {1, 1}, {1, 1}
    };
    private static final double[] CONTRAST = {0, 1};

    @Test void matchesFrozenLimmaVoomEdgeRAndDeseq2Contracts() throws Exception {
        double[][] continuous = matrix("continuous.tsv");
        double[][] counts = matrix("counts.tsv");
        double[][] limma = numeric("limma.tsv", 1, 5);
        double[][] voom = numeric("voom.tsv", 1, 5);
        double[][] edgeR = numeric("edger.tsv", 1, 3);
        double[][] deseq2 = numeric("deseq2.tsv", 1, 3);

        DifferentialFit gaussian = EmpiricalBayesDifferential.fitContinuous(
            continuous, DESIGN, CONTRAST);
        DifferentialFit weighted = EmpiricalBayesDifferential.fitVoom(
            counts, DESIGN, CONTRAST);
        DifferentialFit negativeBinomial =
            EmpiricalBayesDifferential.fitNegativeBinomial(
                counts, DESIGN, CONTRAST);

        double[] javaLimmaEffect = gaussian.results().stream()
            .mapToDouble(DifferentialResult::effect).toArray();
        double[] javaLimmaT = gaussian.results().stream()
            .mapToDouble(DifferentialResult::statistic).toArray();
        double[] javaVoomEffect = weighted.results().stream()
            .mapToDouble(DifferentialResult::effect).toArray();
        double[] javaNbEffect = negativeBinomial.results().stream()
            .mapToDouble(DifferentialResult::log2FoldChange).toArray();
        for (int feature = 0; feature < continuous.length; feature++)
            assertEquals(limma[feature][0], javaLimmaEffect[feature], 2e-12);
        double limmaCorrelation = correlation(javaLimmaT, column(limma, 1));
        double voomCorrelation = correlation(javaVoomEffect, column(voom, 0));
        double edgeCorrelation = correlation(javaNbEffect, column(edgeR, 0));
        double deseqCorrelation = correlation(javaNbEffect, column(deseq2, 0));
        assertTrue(limmaCorrelation > 0.995,
            "moderated t correlation with limma: " + limmaCorrelation);
        assertTrue(voomCorrelation > 0.97,
            "log-fold-change correlation with voom: " + voomCorrelation);
        assertTrue(edgeCorrelation > 0.90,
            "log-fold-change correlation with edgeR: " + edgeCorrelation);
        assertTrue(deseqCorrelation > 0.90,
            "log-fold-change correlation with DESeq2: " + deseqCorrelation);
        for (int feature = 0; feature < 14; feature++) {
            assertEquals(Math.signum(edgeR[feature][0]),
                Math.signum(javaNbEffect[feature]), 0.0);
            assertEquals(Math.signum(deseq2[feature][0]),
                Math.signum(javaNbEffect[feature]), 0.0);
        }
    }

    @Test void fixtureRecordsExactReferenceVersions() throws Exception {
        List<String[]> rows = table("versions.tsv");
        assertTrue(rows.stream().anyMatch(row -> row[0].equals("limma")
            && row[1].equals("3.68.5")));
        assertTrue(rows.stream().anyMatch(row -> row[0].equals("edgeR")
            && row[1].equals("4.10.5")));
        assertTrue(rows.stream().anyMatch(row -> row[0].equals("DESeq2")
            && row[1].equals("1.52.0")));
    }

    private static double[][] matrix(String name) throws IOException {
        List<String[]> rows = table(name);
        double[][] result = new double[rows.size()][rows.get(0).length - 1];
        for (int row = 0; row < rows.size(); row++)
            for (int column = 1; column < rows.get(row).length; column++)
                result[row][column - 1] = Double.parseDouble(rows.get(row)[column]);
        return result;
    }

    private static double[][] numeric(String name, int first, int columns)
            throws IOException {
        List<String[]> rows = table(name);
        double[][] result = new double[rows.size()][columns];
        for (int row = 0; row < rows.size(); row++)
            for (int column = 0; column < columns; column++)
                result[row][column] = Double.parseDouble(rows.get(row)[first + column]);
        return result;
    }

    private static List<String[]> table(String name) throws IOException {
        String resource = "/inference-workflows/" + name;
        var stream = DifferentialReferenceTest.class.getResourceAsStream(resource);
        if (stream == null) throw new IOException("missing test resource: " + resource);
        List<String[]> result = new ArrayList<>();
        try (BufferedReader input = new BufferedReader(new InputStreamReader(
                stream, StandardCharsets.UTF_8))) {
            input.readLine();
            for (String line; (line = input.readLine()) != null;)
                if (!line.isBlank()) result.add(line.split("\\t", -1));
        }
        return result;
    }

    private static double[] column(double[][] values, int column) {
        double[] result = new double[values.length];
        for (int row = 0; row < values.length; row++) result[row] = values[row][column];
        return result;
    }

    private static double correlation(double[] left, double[] right) {
        double leftMean = java.util.Arrays.stream(left).average().orElseThrow();
        double rightMean = java.util.Arrays.stream(right).average().orElseThrow();
        double covariance = 0.0;
        double leftVariance = 0.0;
        double rightVariance = 0.0;
        for (int index = 0; index < left.length; index++) {
            covariance += (left[index] - leftMean) * (right[index] - rightMean);
            leftVariance += (left[index] - leftMean) * (left[index] - leftMean);
            rightVariance += (right[index] - rightMean) * (right[index] - rightMean);
        }
        return covariance / Math.sqrt(leftVariance * rightVariance);
    }
}
