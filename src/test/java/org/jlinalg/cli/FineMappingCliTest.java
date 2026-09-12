/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class FineMappingCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void susieWritesPipsAndColocReadyCredibleEffects() throws Exception {
        Path summary = temporaryDirectory.resolve("summary.tsv");
        Path ld = temporaryDirectory.resolve("ld.tsv");
        Path output = temporaryDirectory.resolve("susie.csv");
        Files.writeString(summary,
            "SNP\tbeta\tse\n"
            + "rs1\t1.2\t0.1\n"
            + "rs2\t0.02\t0.1\n"
            + "rs3\t-0.01\t0.1\n"
            + "rs4\t-1.0\t0.1\n"
            + "rs5\t0.03\t0.1\n");
        Files.writeString(ld,
            "variant_id\trs3\trs1\trs5\trs4\trs2\n"
            + "rs3\t1\t0\t0\t0\t0\n"
            + "rs1\t0\t1\t0\t0\t0\n"
            + "rs5\t0\t0\t1\t0\t0\n"
            + "rs4\t0\t0\t0\t1\t0\n"
            + "rs2\t0\t0\t0\t0\t1\n");

        int status = JLinAlgCli.run(new String[] {
            "susie", "--summary", summary.toString(),
            "--ld", ld.toString(), "--sample-size", "1000",
            "--effects", "2", "--fixed-residual-variance",
            "--tolerance", "1e-8", "--out", output.toString()
        });

        assertEquals(0, status);
        RunLogTest.assertTiming(Path.of(output + ".log"), "complete");
        DelimitedData variants = DelimitedData.read(output);
        assertEquals(5, variants.rows().size());
        assertTrue(variants.header().contains("pip"));
        assertFalse(Files.readString(output).contains("\t"));
        Path effects = Path.of(output + ".effects.tsv");
        DelimitedData credible = DelimitedData.read(effects);
        assertEquals(10, credible.rows().size());
        assertTrue(credible.header().contains("log_bayes_factor"));
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("credible_sets=2"));
    }

    @Test
    void explicitBetaAndSeOverrideAnAutomaticallyRecognizableZColumn()
            throws Exception {
        Path summary = temporaryDirectory.resolve("summary-with-z.tsv");
        Path ld = temporaryDirectory.resolve("ordered-ld.tsv");
        Path output = temporaryDirectory.resolve("explicit-beta.tsv");
        Files.writeString(summary,
            "variant_id\tz\tbeta\tse\n"
            + "rs1\tnot-used\t0.8\t0.1\n"
            + "rs2\tnot-used\t0.0\t0.1\n"
            + "rs3\tnot-used\t-0.7\t0.1\n");
        Files.writeString(ld,
            "1\t0\t0\n0\t1\t0\n0\t0\t1\n");

        int status = JLinAlgCli.run(new String[] {
            "susie", "--summary", summary.toString(),
            "--beta-column", "beta", "--se-column", "se",
            "--ld", ld.toString(), "--sample-size", "500",
            "--effects", "2", "--fixed-residual-variance",
            "--out", output.toString()
        });

        assertEquals(0, status);
        assertEquals(4, Files.readAllLines(output).size());
    }

    @Test
    void colocConsumesSusieEffectTablesAndWritesBothPosteriorLevels()
            throws Exception {
        Path first = effects("trait1", new double[] {12, 0, -1});
        Path second = effects("trait2", new double[] {11, -1, 0});
        Path output = temporaryDirectory.resolve("coloc.csv");

        int status = JLinAlgCli.run(new String[] {
            "coloc", "--trait1", first.toString(),
            "--trait2", second.toString(), "--out", output.toString()
        });

        assertEquals(0, status);
        DelimitedData pairs = DelimitedData.read(output);
        assertEquals(1, pairs.rows().size());
        assertTrue(pairs.header().contains("posterior_h4"));
        assertFalse(Files.readString(output).contains("\t"));
        DelimitedData variants = DelimitedData.read(
            Path.of(output + ".variants.tsv"));
        assertEquals(3, variants.rows().size());
        assertTrue(variants.header().contains(
            "posterior_shared_given_h4"));
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("common_variants=3"));
    }

    @Test
    void colocAcceptsClosedOverlapBoundary() throws Exception {
        Path first = effects("strict-trait1", new double[] {12, 0, -1});
        Path second = effects("strict-trait2", new double[] {11, -1, 0});
        Path output = temporaryDirectory.resolve("strict-coloc.tsv");

        int status = JLinAlgCli.run(new String[] {
            "coloc", "--trait1", first.toString(),
            "--trait2", second.toString(), "--minimum-overlap", "1",
            "--out", output.toString()
        });

        assertEquals(0, status);
        assertTrue(Files.exists(output));
    }

    private Path effects(String name, double[] logBf) throws Exception {
        Path path = temporaryDirectory.resolve(name + ".tsv");
        Files.writeString(path,
            "effect_index\tvariant_id\tlog_bayes_factor\n"
            + "1\trs1\t" + logBf[0] + "\n"
            + "1\trs2\t" + logBf[1] + "\n"
            + "1\trs3\t" + logBf[2] + "\n");
        return path;
    }
}
