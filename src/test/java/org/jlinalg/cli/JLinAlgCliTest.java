/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class JLinAlgCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void ldDatabaseListShowsDownloadChoice() {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream output = new PrintStream(bytes);

        int status = JLinAlgCli.run(new String[] {"ld-db", "list"},
            output, output);

        assertEquals(0, status);
        String text = bytes.toString();
        assertTrue(text.contains("1000g-phase3"));
        assertTrue(text.contains("GRCh37"));
        assertTrue(text.contains("--database 1000g-phase3"));
    }

    @Test
    void ldDatabaseDownloadWithoutChoiceReturnsActionableError() {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream output = new PrintStream(bytes);

        int status = JLinAlgCli.run(new String[] {"ld-db", "download"},
            output, output);

        assertEquals(2, status);
        String text = bytes.toString();
        assertTrue(text.contains("no LD database was specified"));
        assertTrue(text.contains("--database is required"));
        assertTrue(text.contains("Database choices:"));
        assertTrue(text.contains("ld-db download"));
    }

    @Test
    void blankLdDatabaseLocationUsesCurrentDirectory() {
        AtomicReference<Path> installedAt = new AtomicReference<>();
        PrintStream output = new PrintStream(OutputStream.nullOutputStream());

        int status = LdDatabaseCli.run(new String[] {
            "download", "--database", "1000g-phase3", "--location", ""
        }, output, output, temporaryDirectory,
            (database, location, progress) -> installedAt.set(location));

        assertEquals(0, status);
        assertEquals(temporaryDirectory.toAbsolutePath().normalize(),
            installedAt.get());
    }

    @Test
    void phenotypeOnlyOlsWritesAllTermsFdrLogAndManifest() throws Exception {
        Path phenotype = temporaryDirectory.resolve("phenotype.tsv");
        Path output = temporaryDirectory.resolve("fit.tsv");
        Files.writeString(phenotype,
            "IID\ty\tage\n"
            + "S1\t1.0\t20\n"
            + "S2\t2.0\t30\n"
            + "S3\t3.2\t40\n"
            + "S4\t3.8\t50\n"
            + "S5\t5.1\t60\n");

        int status = JLinAlgCli.run(new String[] {
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "y ~ age", "--out", output.toString()
        });

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(output);
        assertEquals(3, lines.size());
        assertTrue(lines.get(0).endsWith("fdr_bh"));
        assertTrue(lines.stream().anyMatch(line ->
            line.startsWith("ok\tage\t")));
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("resolved_model=ols"));
        assertTrue(Files.readString(Path.of(output + ".manifest.json"))
            .contains("\"model\": \"ols\""));
    }

    @Test
    void streamsEwasRowsWithTransformAnnotationAndBh() throws Exception {
        Path phenotype = temporaryDirectory.resolve("phenotype.csv");
        Path omics = temporaryDirectory.resolve("methylation.csv");
        Path annotation = temporaryDirectory.resolve("annotation.tsv");
        Path output = temporaryDirectory.resolve("ewas.tsv");
        Files.writeString(phenotype,
            "IID,trait,age\n"
            + "S3,3.2,40\nS1,1.0,20\nS4,3.8,50\n"
            + "S2,2.0,30\nS6,6.0,70\nS5,5.1,60\n");
        Files.writeString(omics,
            "probe,S1,S2,S3,S4,S5,S6\n"
            + "cg12345678,0.10,0.20,0.30,0.40,0.50,0.60\n"
            + "cg87654321,0.60,0.50,0.40,0.30,0.20,0.10\n");
        Files.writeString(annotation,
            "probe\tchr\tstart\tgene\n"
            + "cg12345678\t1\t100\tGENE1\n"
            + "cg87654321\t2\t200\tGENE2\n");

        int status = JLinAlgCli.run(new String[] {
            "--omics", omics.toString(), "--pheno", phenotype.toString(),
            "--id", "IID", "--formula", "trait ~ age + <omics>",
            "--transform", "<omics>=mvalue(epsilon=1e-6)|zscore()",
            "--annot", annotation.toString(), "--annot-id", "probe",
            "--annot-cols", "chr,start,gene", "--block-size", "1",
            "--threads", "1", "--out", output.toString()
        });

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(output);
        assertEquals(3, lines.size());
        assertTrue(lines.get(0).contains("annot_gene"));
        assertTrue(lines.get(0).endsWith("fdr_bh"));
        assertTrue(lines.get(1).contains("\tewas\tcg12345678\t"));
        assertTrue(lines.get(1).contains("\tGENE1\t"));
        assertFalse(Files.exists(Path.of(output + ".partial")));
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("block_size=1"));
    }

    @Test
    void gwasOutputIsAltOrientedAndIncludesStratifiedHwe() throws Exception {
        Path phenotype = temporaryDirectory.resolve("case-control.tsv");
        Path variants = temporaryDirectory.resolve("variants.tsv");
        Path output = temporaryDirectory.resolve("gwas.tsv");
        Files.writeString(phenotype,
            "IID\tcase\tage\n"
            + "S1\tNo\t20\nS2\tNo\t30\nS3\tNo\t40\n"
            + "S4\tYes\t50\nS5\tYes\t60\nS6\tYes\t70\n");
        Files.writeString(variants,
            "id\tchr\tposition\tref\talt\tS1\tS2\tS3\tS4\tS5\tS6\n"
            + "rs123\t1\t100\tA\tG\t0\t0\t1\t1\t2\t2\n"
            + "rs456\t1\t200\tC\tT\t0\t1\t0\t1\t0\t1\n");

        int status = JLinAlgCli.run(new String[] {
            "--omics", variants.toString(), "--pheno", phenotype.toString(),
            "--id", "IID", "--formula", "case ~ <omics>",
            "--family", "binomial", "--threads", "1", "--block-size", "1",
            "--out", output.toString()
        });

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(output);
        assertEquals(3, lines.size());
        String[] header = lines.get(0).split("\t", -1);
        String[] row = lines.get(1).split("\t", -1);
        assertEquals("G", row[index(header, "effect_allele")]);
        assertFalse(row[index(header, "hwe_p_all")].isEmpty());
        assertFalse(row[index(header, "hwe_p_cases")].isEmpty());
        assertFalse(row[index(header, "hwe_p_controls")].isEmpty());
        assertEquals("t_approx", row[index(header, "statistic_type")]);
    }

    @Test
    void grmSelectsRemlAndAlignsThroughIndividualId() throws Exception {
        Path phenotype = temporaryDirectory.resolve("repeated.tsv");
        Path grm = temporaryDirectory.resolve("grm.tsv");
        Path output = temporaryDirectory.resolve("reml.tsv");
        Files.writeString(phenotype,
            "observation\tIID\ty\tage\n"
            + "O3\tS3\t4.0\t30\nO1\tS1\t2.0\t10\n"
            + "O2\tS2\t3.1\t20\nO6\tS6\t6.1\t60\n"
            + "O4\tS4\t4.7\t40\nO5\tS5\t5.5\t50\n"
            + "O9\tS9\t8.8\t90\nO7\tS7\t7.0\t70\n"
            + "O8\tS8\t7.9\t80\nO12\tS12\t11.7\t120\n"
            + "O10\tS10\t10.2\t100\nO11\tS11\t10.8\t110\n");
        writeBlockGrm(grm, 12);

        int status = JLinAlgCli.run(new String[] {
            "--pheno", phenotype.toString(), "--id", "observation",
            "--individual-id", "IID", "--formula", "y ~ age",
            "--grm", grm.toString(), "--out", output.toString()
        });

        assertEquals(0, status);
        assertTrue(Files.exists(output));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("resolved_model=lmm"));
        assertTrue(log.contains("grm_match_column=IID"));
        String manifest = Files.readString(Path.of(output + ".manifest.json"));
        assertTrue(manifest.contains("\"grm\":"));
    }

    @Test
    void gwasGrmUsesNullModelWithoutFormulaRandomTerm() throws Exception {
        Path phenotype = temporaryDirectory.resolve("quantitative.tsv");
        Path variants = temporaryDirectory.resolve("grm-variants.tsv");
        Path grm = temporaryDirectory.resolve("scan-grm.tsv");
        Path output = temporaryDirectory.resolve("grm-gwas.tsv");
        StringBuilder observations = new StringBuilder("IID\ty\tage\n");
        for (int id = 1; id <= 12; id++)
            observations.append('S').append(id).append('\t')
                .append(1.5 + id * 0.3 + (id % 3) * 0.4).append('\t')
                .append(20 + id).append('\n');
        Files.writeString(phenotype, observations);
        StringBuilder genotypes = new StringBuilder(
            "id\tchr\tposition\tref\talt");
        for (int id = 1; id <= 12; id++)
            genotypes.append("\tS").append(id);
        genotypes.append('\n').append("rs101\t1\t101\tA\tG");
        for (int id = 1; id <= 12; id++)
            genotypes.append('\t').append(id % 3);
        genotypes.append('\n').append("rs202\t2\t202\tC\tT");
        for (int id = 1; id <= 12; id++)
            genotypes.append('\t').append((id / 2) % 3);
        genotypes.append('\n');
        Files.writeString(variants, genotypes);
        writeBlockGrm(grm, 12);

        int status = JLinAlgCli.run(new String[] {
            "--omics", variants.toString(),
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "y ~ age + <omics>", "--grm", grm.toString(),
            "--threads", "1", "--block-size", "1",
            "--out", output.toString()
        });

        assertEquals(0, status);
        assertEquals(3, Files.readAllLines(output).size());
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("resolved_model=lmm"));
    }

    private static void writeBlockGrm(Path path, int size) throws Exception {
        StringBuilder matrix = new StringBuilder("IID");
        for (int id = 1; id <= size; id++) matrix.append("\tS").append(id);
        matrix.append('\n');
        for (int row = 1; row <= size; row++) {
            matrix.append('S').append(row);
            for (int column = 1; column <= size; column++) {
                double value = row == column ? 1.0
                    : (row - 1) / 3 == (column - 1) / 3 ? 0.35 : 0.0;
                matrix.append('\t').append(value);
            }
            matrix.append('\n');
        }
        Files.writeString(path, matrix);
    }

    private static int index(String[] values, String target) {
        for (int index = 0; index < values.length; index++)
            if (values[index].equals(target)) return index;
        throw new AssertionError("column not found: " + target);
    }
}
