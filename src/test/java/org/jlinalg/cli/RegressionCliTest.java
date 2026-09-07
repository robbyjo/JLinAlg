/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;

class RegressionCliTest {
    @Test
    void betaRegressionCommandFitsNumericCsv() throws Exception {
        Path input = Files.createTempFile("jlinalg-beta-", ".csv");
        Files.writeString(input, "y,x\n0.20,0\n0.30,1\n0.42,2\n0.55,3\n0.68,4\n0.78,5\n");
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        int status = JLinAlgCli.run(new String[] {"beta-regression", "--input",
            input.toString(), "--response", "y", "--mean", "x"},
            new PrintStream(output, true, StandardCharsets.UTF_8),
            new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
        assertEquals(0, status);
        assertTrue(output.toString(StandardCharsets.UTF_8).contains("(Intercept)"));
        assertTrue(output.toString(StandardCharsets.UTF_8).contains("converged"));
        Files.deleteIfExists(input);
    }

    @Test
    void penalizedRegressionCommandSupportsElasticNetAndCv() throws Exception {
        Path input = Files.createTempFile("jlinalg-penalized-", ".tsv");
        Files.writeString(input, "y\tx1\tx2\n1\t0\t1\n2\t1\t0\n3\t2\t1\n4\t3\t0\n5\t4\t1\n6\t5\t0\n");
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        int status = JLinAlgCli.run(new String[] {"penalized-regression",
            "--input", input.toString(), "--response", "y", "--predictors", "x1,x2",
            "--model", "elastic-net", "--alpha", "0.5", "--lambda-grid", "0.1,0.01",
            "--cv-folds", "2"},
            new PrintStream(output, true, StandardCharsets.UTF_8),
            new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
        assertEquals(0, status);
        assertTrue(output.toString(StandardCharsets.UTF_8).contains("elastic-net"));
        assertTrue(output.toString(StandardCharsets.UTF_8).contains("x1"));
        Files.deleteIfExists(input);
    }
}
