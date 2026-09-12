/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class GrmCliTest {
    @TempDir Path directory;

    private int run(String... args) {
        return JLinAlgCli.run(args, new PrintStream(new ByteArrayOutputStream()),
            new PrintStream(new ByteArrayOutputStream()));
    }

    @Test void analyticMissingDosagesAndUnequalBlocksRoundTrip() throws Exception {
        Path input = directory.resolve("dosage.tsv");
        Files.writeString(input, "id\tsample_id\tb\tc\td\n"
            + "v1\t0\t1\t2\tNA\nv2\t0\t0\t1\t1\nmono\t0\t0\t0\t0\n"
            + "lowcall\tNA\tNA\tNA\t1\nv5\t0\t0\t0\t2\n");
        // Independently derived sum of three standardized outer products / 3.
        double[] expected = {10./9,4./9,-2./3,-8./9, 4./9,4./9,0,-8./9,
            -2./3,0,10./9,-4./9, -8./9,-8./9,-4./9,20./9};
        for (int blockSize : new int[]{1, 2, 5}) {
            Path output = directory.resolve("grm" + blockSize + (blockSize == 2 ? ".csv" : ".tsv"));
            assertEquals(0, run("grm", "--genotypes", input.toString(), "--out", output.toString(),
                "--maf", "0", "--call-rate", "0.5", "--block-size", Integer.toString(blockSize)));
            var loaded = GrmReader.read(output).matrix();
            assertEquals(List.of("sample_id", "b", "c", "d"), loaded.sampleIds());
            assertArrayEquals(expected, loaded.relationshipMatrix(), 1e-13);
            String log = Files.readString(Path.of(output + ".log"));
            for (String field : List.of("variants_considered=5", "variants_used=3", "variants_excluded=2",
                    "started=", "finished=", "elapsed=", "status=complete")) assertTrue(log.contains(field), field);
            String original = Files.readString(output);
            assertEquals(2, run("grm", "--genotypes", input.toString(), "--out", output.toString()));
            assertEquals(original, Files.readString(output));
        }
    }

    @Test void emptyAndInvalidInputsDoNotPublishMatrices() throws Exception {
        Path input = directory.resolve("bad.tsv"), output = directory.resolve("failed.tsv");
        Files.writeString(input, "id\ta\tb\nv\t0\t0\n");
        assertEquals(2, run("grm", "--genotypes", input.toString(), "--out", output.toString()));
        assertFalse(Files.exists(output));
        assertTrue(Files.readString(Path.of(output + ".log")).contains("status=failed"));
        Files.writeString(input, "id\ta\tb\nv\t0\t3\n");
        assertEquals(2, run("grm", "--genotypes", input.toString(), "--out", directory.resolve("invalid.tsv").toString()));
        assertFalse(Files.exists(directory.resolve("invalid.tsv")));
        assertEquals(2, run("grm", "--genotypes", input.toString(), "--out", directory.resolve("options.tsv").toString(), "--block-size", "0"));
        assertFalse(Files.exists(directory.resolve("options.tsv.log")));
        assertEquals(0, run("grm", "--help"));
    }

    @Test void vcfConstructionFeedsTheMixedModelCli() throws Exception {
        Path output = directory.resolve("cohort.tsv");
        assertEquals(0, run("grm", "--genotypes", "examples/rare-meta/cohort-a.vcf", "--out", output.toString()));
        assertEquals(6, GrmReader.read(output).matrix().samples());
        Path fit = directory.resolve("model.tsv");
        assertEquals(0, run("--pheno", "examples/rare-meta/cohort-a.tsv", "--id", "sample",
            "--formula", "trait ~ 1", "--grm", output.toString(), "--backend", "cpu", "--out", fit.toString()));
        assertTrue(Files.size(fit) > 0);
    }
}
