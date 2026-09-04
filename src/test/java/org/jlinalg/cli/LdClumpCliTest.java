/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class LdClumpCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void clumpsCanonicalInstrumentFileWithDownloadedDatabaseLayout()
            throws Exception {
        Path database = temporaryDirectory.resolve("ld");
        Path prefix = database.resolve("panels/EUR/genotypes");
        Files.createDirectories(prefix.getParent());
        PlinkTestReference.write(prefix,
            new PlinkTestReference.Variant("1", "rsLead", 10_000,
                new int[] {0, 0, 1, 1, 2, 2, 0, 2}),
            new PlinkTestReference.Variant("1", "rsLinked", 20_000,
                new int[] {0, 0, 1, 1, 2, 2, 0, 2}),
            new PlinkTestReference.Variant("2", "rsOther", 20_000,
                new int[] {0, 1, 2, 0, 1, 2, 1, 0}));
        Files.writeString(database.resolve("jlinalg-ld-reference.json"),
            "{\"format\":\"jlinalg-ld-reference\",\"format_version\":1,"
                + "\"panels\":[{\"id\":\"EUR\"}]}\n");
        Path instruments = temporaryDirectory.resolve("instruments.tsv");
        Path clumped = temporaryDirectory.resolve("clumped.tsv");
        Files.writeString(instruments,
            "Phenotype\tSNP\tbeta\tpval\n"
                + "BMI\trsLinked\t0.2\t2e-8\n"
                + "BMI\trsLead\t0.1\t1e-9\n"
                + "BMI\trsOther\t0.3\t3e-8\n"
                + "BMI\trsAbsent\t0.4\t4e-8\n");
        ByteArrayOutputStream messages = new ByteArrayOutputStream();
        PrintStream output = new PrintStream(messages);

        int status = LdClumpCli.run(new String[] {
            "--database", database.toString(),
            "--instrument", instruments.toString(),
            "--ld-threshold", "0.001", "--output", clumped.toString()
        }, output, output, temporaryDirectory);

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(clumped);
        assertEquals(List.of(
            "Phenotype\tSNP\tbeta\tpval",
            "BMI\trsLead\t0.1\t1e-9",
            "BMI\trsOther\t0.3\t3e-8"), lines);
        String report = messages.toString();
        assertTrue(report.contains("removed 1 for LD, 1 absent"));
        assertTrue(report.contains("using EUR"));
    }

    @Test
    void rejectsDirectoryWithoutInstalledManifest() {
        ByteArrayOutputStream messages = new ByteArrayOutputStream();
        PrintStream output = new PrintStream(messages);

        int status = LdClumpCli.run(new String[] {
            "--database", temporaryDirectory.toString(),
            "--instrument", "missing.tsv", "--output", "out.tsv"
        }, output, output, temporaryDirectory);

        assertEquals(1, status);
        assertTrue(messages.toString().contains("manifest is absent"));
    }
}
