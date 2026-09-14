/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class ProbitCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void fitsBinaryProbitThroughTheGeneralCli() throws Exception {
        Path phenotype = temporaryDirectory.resolve("probit.csv");
        Path output = temporaryDirectory.resolve("probit-results.csv");
        Files.writeString(phenotype,
            "IID,case,x,z\n"
            + "S1,0,-2,0\nS2,0,-1.7,1\nS3,0,-1.4,0\nS4,0,-1.1,1\n"
            + "S5,1,-.8,0\nS6,0,-.5,1\nS7,0,-.2,0\nS8,1,.1,1\n"
            + "S9,0,.4,0\nS10,1,.7,1\nS11,1,1,0\nS12,0,1.3,1\n"
            + "S13,1,1.6,0\nS14,1,1.9,1\nS15,1,2.2,0\nS16,1,2.5,1\n");

        int status = JLinAlgCli.run(new String[] {
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "case ~ x + z", "--family", "probit",
            "--link", "probit", "--out", output.toString()
        });

        assertEquals(0, status);
        List<String> rows = Files.readAllLines(output);
        assertTrue(rows.size() >= 2);
        String manifest = Files.readString(Path.of(output + ".manifest.json"));
        assertTrue(manifest.contains("\"family\""));
        assertTrue(manifest.contains("\"probit\""));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("statistic_type=z"));
        assertTrue(log.contains("df_method=asymptotic"));
    }

    @Test
    void estimatedDispersionMetadataDistinguishesFitFromScan()
            throws Exception {
        Path phenotype = temporaryDirectory.resolve("gaussian.csv");
        Path variants = temporaryDirectory.resolve("variants.tsv");
        Path fitOutput = temporaryDirectory.resolve("gaussian-fit.csv");
        Path scanOutput = temporaryDirectory.resolve("gaussian-scan.tsv");
        Files.writeString(phenotype,
            "IID,y,x\n"
            + "S1,1.0,0\nS2,2.1,1\nS3,1.8,2\nS4,3.7,3\n"
            + "S5,4.2,4\nS6,5.8,5\nS7,5.5,6\nS8,7.1,7\n");
        Files.writeString(variants,
            "id\tchr\tposition\tref\talt\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\n"
            + "rs1\t1\t100\tA\tG\t0\t1\t2\t0\t1\t2\t0\t1\n");

        assertEquals(0, JLinAlgCli.run(new String[] {
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "y ~ x", "--model", "glm",
            "--family", "gaussian", "--out", fitOutput.toString()
        }));
        String fitLog = Files.readString(Path.of(fitOutput + ".log"));
        assertTrue(fitLog.contains("statistic_type=t"));
        assertTrue(fitLog.contains("df_method=residual"));

        assertEquals(0, JLinAlgCli.run(new String[] {
            "--omics", variants.toString(),
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "y ~ x + <omics>", "--model", "glm",
            "--family", "gaussian", "--threads", "1",
            "--block-size", "1", "--out", scanOutput.toString()
        }));
        String scanLog = Files.readString(Path.of(scanOutput + ".log"));
        assertTrue(scanLog.contains("statistic_type=t_approx"));
        assertTrue(scanLog.contains("df_method=residual-approximation"));
    }
}
