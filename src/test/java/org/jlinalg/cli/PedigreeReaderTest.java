/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import java.nio.file.Files;
import java.nio.file.Path;
import org.jlinalg.pedigree.Pedigree;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PedigreeReaderTest {
    @TempDir Path temporaryDirectory;

    @Test
    void sparseReaderInbreedingMatchesTabularRelationship() throws Exception {
        Path input = temporaryDirectory.resolve("pedigree.tsv");
        Files.writeString(input,
            "id\tsire\tdam\n"
            + "E\tC\tD\n"
            + "A\t0\t0\n"
            + "D\tA\tB\n"
            + "B\t0\t0\n"
            + "C\tA\tB\n");

        PedigreeReader.Loaded loaded =
            PedigreeReader.read(input, "id", "sire", "dam", null);
        double[] expected = Pedigree.of(loaded.individuals())
            .inbreedingCoefficients();

        assertArrayEquals(expected, loaded.inbreedingCoefficients(), 1e-12);
    }
}
