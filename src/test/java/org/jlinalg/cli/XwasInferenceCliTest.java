/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

final class XwasInferenceCliTest {
    @TempDir Path temporaryDirectory;

    @Test void newCommandsAreDispatchedAndAdvertiseHelp() {
        for (String command : new String[] {"differential", "ewas-regions",
                "multiple-test", "multiple-impute", "mi-pool"}) {
            ByteArrayOutputStream bytes = new ByteArrayOutputStream();
            assertEquals(0, JLinAlgCli.run(new String[] {command, "--help"},
                new PrintStream(bytes), System.err));
            assertTrue(bytes.toString().contains("Usage: jlinalg " + command));
        }
    }

    @Test void differentialAndRegionCommandsPublishAuditableTables() throws Exception {
        Path omics = temporaryDirectory.resolve("omics.tsv");
        Path phenotype = temporaryDirectory.resolve("phenotype.tsv");
        Files.writeString(omics, "feature\ts1\ts2\ts3\ts4\ts5\ts6\n"
            + "g1\t1\t1.2\t0.8\t3\t3.1\t2.9\n"
            + "g2\t2\t2.1\t1.9\t2.1\t2\t2.2\n"
            + "g3\t4\t4.3\t3.7\t5\t5.2\t4.8\n"
            + "g4\t8\t7\t9\t8.1\t7.2\t8.9\n");
        Files.writeString(phenotype, "sample\tgroup\tage\n"
            + "s1\tcontrol\t40\n" + "s2\tcontrol\t41\n"
            + "s3\tcontrol\t39\n" + "s4\tcase\t40\n"
            + "s5\tcase\t41\n" + "s6\tcase\t39\n");
        Path differential = temporaryDirectory.resolve("differential.tsv");
        assertEquals(0, JLinAlgCli.run(new String[] {"differential",
            "--method", "limma", "--omics", omics.toString(), "--pheno",
            phenotype.toString(), "--id", "sample", "--group", "group",
            "--out", differential.toString()}));
        assertTrue(Files.readString(differential).startsWith(
            "feature_id\teffect\tlog2_fold_change"));
        assertTrue(Files.readString(temporaryDirectory.resolve(
            "differential.metadata.tsv")).contains("aligned_samples\t6"));

        Path probes = temporaryDirectory.resolve("probes.tsv");
        Files.writeString(probes,
            "probe_id\tchromosome\tposition\teffect\tp\n"
            + "a\t1\t100\t0.2\t0.01\n" + "b\t1\t150\t0.3\t0.02\n"
            + "c\t1\t220\t0.1\t0.05\n" + "d\t2\t1000\t-0.2\t0.1\n"
            + "e\t2\t1060\t-0.3\t0.03\n" + "f\t2\t1120\t-0.1\t0.04\n");
        Path regions = temporaryDirectory.resolve("regions.tsv");
        assertEquals(0, JLinAlgCli.run(new String[] {"ewas-regions",
            "--input", probes.toString(), "--genome-build", "GRCh38",
            "--out", regions.toString()}));
        String regionText = Files.readString(regions);
        assertTrue(regionText.contains("probes_per_kb"));
        assertTrue(regionText.contains("GRCh38\tchr1\t100\t220\t3"));
    }

    @Test void familyAccountingImputationAndRubinPoolingAreRunnable() throws Exception {
        Path family = temporaryDirectory.resolve("family.tsv");
        Files.writeString(family, "id\tparent\tp\tstatus\n"
            + "gene\t\t0.001\tok\n" + "tissue\tgene\t0.01\tok\n"
            + "trait\ttissue\tNA\tfailed\n");
        Path adjusted = temporaryDirectory.resolve("adjusted.tsv");
        assertEquals(0, JLinAlgCli.run(new String[] {"multiple-test",
            "--method", "hierarchy", "--input", family.toString(),
            "--parent", "parent", "--status", "status", "--out",
            adjusted.toString()}));
        assertTrue(Files.readString(adjusted).contains("failed\t1.0"));

        Path incomplete = temporaryDirectory.resolve("incomplete.tsv");
        Files.writeString(incomplete, "sample\tage\tcase\tsite\n"
            + "s1\t20\t0\tA\n" + "s2\t21\t0\tA\n"
            + "s3\t22\t0\tB\n" + "s4\tNA\t1\tB\n"
            + "s5\t30\t1\tNA\n" + "s6\t31\t1\tB\n"
            + "s7\t32\tNA\tA\n" + "s8\t33\t1\tB\n");
        Path prefix = temporaryDirectory.resolve("completed");
        assertEquals(0, JLinAlgCli.run(new String[] {"multiple-impute",
            "--input", incomplete.toString(), "--id", "sample", "--types",
            "age:continuous,case:binary,site:categorical", "--imputations", "2",
            "--iterations", "3", "--seed", "7", "--out", prefix.toString()}));
        assertTrue(Files.exists(temporaryDirectory.resolve("completed.imp1.tsv")));
        assertTrue(Files.readString(temporaryDirectory.resolve(
            "completed.diagnostics.tsv")).contains("between_chain_variance"));

        Path estimates = temporaryDirectory.resolve("estimates.tsv");
        Files.writeString(estimates,
            "parameter\timputation\testimate\tvariance\n"
            + "group\t1\t1\t4\n" + "group\t2\t2\t4\n"
            + "group\t3\t3\t4\n");
        Path pooled = temporaryDirectory.resolve("pooled.tsv");
        assertEquals(0, JLinAlgCli.run(new String[] {"mi-pool",
            "--input", estimates.toString(), "--complete-df", "20",
            "--out", pooled.toString()}));
        assertTrue(Files.readString(pooled).contains("group\t3\t2.0\t4.0\t1.0"));
    }
}
