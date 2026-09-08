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

final class MrEstimatorExtensionTest {
    @Test
    void generalizedAndConditionalMethodsAcceptLdMatrix() throws Exception {
        Path input = Files.createTempFile("jlinalg-mr-", ".tsv"); Path ld = Files.createTempFile("jlinalg-ld-", ".tsv");
        Files.writeString(input, "variant_id\tbeta_exposure\tse_exposure\tbeta_outcome\tse_outcome\nrs1\t0.2\t0.05\t0.1\t0.04\nrs2\t0.3\t0.06\t0.2\t0.05\nrs3\t0.4\t0.07\t0.25\t0.06\nrs4\t0.5\t0.08\t0.3\t0.07\n");
        Files.writeString(ld, "1\t0.1\t0\t0\n0.1\t1\t0.1\t0\n0\t0.1\t1\t0.1\n0\t0\t0.1\t1\n");
        ByteArrayOutputStream output = new ByteArrayOutputStream(); int status = JLinAlgCli.run(new String[] {"mr-estimate", "--input", input.toString(), "--method", "conditional", "--ld", ld.toString()}, new PrintStream(output, true, StandardCharsets.UTF_8), new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
        assertEquals(0, status); assertTrue(output.toString(StandardCharsets.UTF_8).contains("IVW_GENERALIZED_FIXED")); Files.deleteIfExists(input); Files.deleteIfExists(ld);
    }
}
