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
        Path output = temporaryDirectory.resolve("ewas.csv");
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
        DelimitedData result = DelimitedData.read(output);
        assertEquals(2, result.rows().size());
        assertTrue(result.header().contains("annot_gene"));
        assertEquals("fdr_bh",
            result.header().get(result.header().size() - 1));
        String[] first = result.rows().get(0);
        assertEquals("ok", first[result.column("status")]);
        assertEquals("cg12345678", first[result.column("id")]);
        assertEquals("GENE1", first[result.column("annot_gene")]);
        assertFalse(Files.readString(output).contains("\t"));
        assertFalse(Files.exists(Path.of(output + ".partial")));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("block_size=1"));
        assertTrue(log.contains("association_chunk_size=1"));
        assertTrue(log.contains("scan_worker_capacity=1"));
        assertTrue(log.contains("omics_type=ewas"));
        assertTrue(log.contains("output_format=csv"));
        String manifest = Files.readString(
            Path.of(output + ".manifest.json"));
        assertTrue(manifest.contains(
            "\"association_chunk_size\": \"1\""));
        assertTrue(manifest.contains(
            "\"scan_worker_capacity\": \"1\""));
    }

    @Test
    void intersectsMismatchedOmicsAndPhenotypeSamplesAndReportsCounts()
            throws Exception {
        Path phenotype = temporaryDirectory.resolve("mismatched-phenotype.csv");
        Path omics = temporaryDirectory.resolve("mismatched-expression.csv");
        Path output = temporaryDirectory.resolve("aligned.tsv");
        Files.writeString(phenotype,
            "IID,trait,age\n"
            + "S3,3.2,40\nS1,1.0,20\nS4,,50\n"
            + "S2,2.0,30\nS7,7.0,80\nS5,5.1,60\n");
        Files.writeString(omics,
            "gene,S1,S2,S3,S4,S5,S6\n"
            + "GENE1,0.10,0.22,0.31,0.43,0.54,0.65\n"
            + "GENE2,0.61,0.52,0.44,0.33,0.21,0.12\n");
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream console = new PrintStream(bytes);

        int status = JLinAlgCli.run(new String[] {
            "--omics", omics.toString(), "--pheno", phenotype.toString(),
            "--id", "IID", "--formula", "trait ~ age + <omics>",
            "--threads", "1", "--block-size", "1",
            "--out", output.toString()
        }, console, console);

        assertEquals(0, status);
        assertEquals(3, Files.readAllLines(output).size());
        assertTrue(bytes.toString().contains(
            "Aligned samples: 5 (omics=6, phenotype=6, omics-only=1, "
                + "phenotype-only=1)"));
        assertTrue(bytes.toString().contains(
            "Analysis samples: 4 (phenotype-missing omitted=1)"));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("omics_samples=6"));
        assertTrue(log.contains("phenotype_samples=6"));
        assertTrue(log.contains("aligned_samples=5"));
        assertTrue(log.contains("omics_only_samples=1"));
        assertTrue(log.contains("phenotype_only_samples=1"));
        assertTrue(log.contains("phenotype_missing_samples_omitted=1"));
        assertTrue(log.contains("analysis_samples=4"));
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
        assertFalse(List.of(header).contains("statistic_type"));
        assertTrue(Files.readString(Path.of(output + ".log"))
            .contains("statistic_type=t_approx"));
    }

    @Test
    void maximumMafSelectsRareVariantsAndReportsFilteredCommonVariants()
            throws Exception {
        Path phenotype = temporaryDirectory.resolve("rare-pheno.tsv");
        Path variants = temporaryDirectory.resolve("rare-variants.tsv");
        Path output = temporaryDirectory.resolve("rare-scan.tsv");
        Files.writeString(phenotype,
            "IID\ty\tage\n"
            + "S1\t1.0\t20\nS2\t2.2\t25\nS3\t1.7\t30\n"
            + "S4\t3.5\t35\nS5\t2.8\t40\nS6\t4.4\t45\n"
            + "S7\t3.9\t50\nS8\t5.3\t55\nS9\t4.8\t60\n"
            + "S10\t6.1\t65\n");
        Files.writeString(variants,
            "id\tchr\tposition\tref\talt\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10\n"
            + "rs100\t1\t100\tA\tG\t0\t0\t0\t1\t0\t0\t0\t0\t0\t0\n"
            + "rs200\t1\t200\tC\tT\t0\t1\t2\t0\t1\t2\t0\t1\t2\t1\n");

        int status = JLinAlgCli.run(new String[] {
            "--omics", variants.toString(), "--pheno", phenotype.toString(),
            "--id", "IID", "--formula", "y ~ age + <omics>",
            "--max-maf", "0.1", "--max-mac", "2",
            "--threads", "1", "--block-size", "1",
            "--out", output.toString()
        });

        assertEquals(0, status);
        DelimitedData result = DelimitedData.read(output);
        String[] rare = result.rows().stream()
            .filter(row -> row[result.column("id")].equals("rs100"))
            .findFirst().orElseThrow();
        String[] common = result.rows().stream()
            .filter(row -> row[result.column("id")].equals("rs200"))
            .findFirst().orElseThrow();
        assertEquals("ok", rare[result.column("status")]);
        assertEquals("filtered", common[result.column("status")]);
        assertTrue(common[result.column("filter_reason")]
            .contains("ABOVE_MAXIMUM_MAF"));
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

    @Test
    void numericMixedScanDefaultsToPerFeatureRemlRefits() throws Exception {
        Path phenotype = temporaryDirectory.resolve("mixed-expression-pheno.tsv");
        Path omics = temporaryDirectory.resolve("mixed-expression.tsv");
        Path output = temporaryDirectory.resolve("mixed-expression-results.tsv");
        StringBuilder observations = new StringBuilder("IID\ty\tx\tbatch\n");
        StringBuilder features = new StringBuilder("gene");
        for (int row = 0; row < 48; row++) {
            int group = row / 4;
            double x = (row % 7) - 3;
            double gene = ((row * 5) % 13 - 6) / 3.0;
            double noise = ((row * 17) % 9 - 4) / 20.0;
            observations.append('S').append(row).append('\t')
                .append(2.0 + 0.25 * x + 0.7 * gene
                    + (group % 4 - 1.5) * 0.8 + noise)
                .append('\t').append(x).append("\tB").append(group)
                .append('\n');
            features.append("\tS").append(row);
        }
        features.append('\n').append("GENE1");
        for (int row = 0; row < 48; row++)
            features.append('\t').append(((row * 5) % 13 - 6) / 3.0);
        features.append('\n').append("GENE2");
        for (int row = 0; row < 48; row++)
            features.append('\t').append(((row * 7) % 17 - 8) / 4.0);
        features.append('\n');
        Files.writeString(phenotype, observations);
        Files.writeString(omics, features);

        int status = JLinAlgCli.run(new String[] {
            "--omics", omics.toString(), "--omics-type", "expression",
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "y ~ x + <omics> + (1|batch)",
            "--threads", "2", "--block-size", "2",
            "--backend", "cpu", "--out", output.toString()
        });

        assertEquals(0, status);
        assertEquals(3, Files.readAllLines(output).size());
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("resolved_model=lmm"));
        assertTrue(log.contains("variance_components=refit"));
        assertTrue(log.contains("mixed_fit=exact-reml-refit"));
        assertFalse(log.contains("P3D"));
    }

    @Test
    void numericNonGaussianMixedScanUsesLaplaceGlmm() throws Exception {
        Path phenotype = temporaryDirectory.resolve("mixed-binary-pheno.tsv");
        Path omics = temporaryDirectory.resolve("mixed-binary-expression.tsv");
        Path output = temporaryDirectory.resolve("mixed-binary-results.tsv");
        StringBuilder observations =
            new StringBuilder("IID\tcase\tx\tbatch\n");
        StringBuilder features = new StringBuilder("gene");
        for (int row = 0; row < 60; row++) {
            int group = row / 5;
            double x = (row % 9) - 4;
            int outcome = ((row * 11 + group * 3) % 17)
                < 7 + (x > 0 ? 2 : 0) ? 1 : 0;
            observations.append('S').append(row).append('\t')
                .append(outcome).append('\t').append(x)
                .append("\tB").append(group).append('\n');
            features.append("\tS").append(row);
        }
        features.append('\n').append("GENE1");
        for (int row = 0; row < 60; row++)
            features.append('\t').append(((row * 5) % 19 - 9) / 4.0);
        features.append('\n');
        Files.writeString(phenotype, observations);
        Files.writeString(omics, features);

        int status = JLinAlgCli.run(new String[] {
            "--omics", omics.toString(), "--omics-type", "expression",
            "--pheno", phenotype.toString(), "--id", "IID",
            "--formula", "case ~ x + <omics> + (1|batch)",
            "--family", "binomial", "--threads", "1",
            "--block-size", "1", "--backend", "cpu",
            "--out", output.toString()
        });

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(output);
        assertEquals(2, lines.size());
        String[] header = lines.get(0).split("\t", -1);
        assertFalse(List.of(header).contains("statistic_type"));
        assertFalse(List.of(header).contains("df_method"));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("resolved_model=glmm"));
        assertTrue(log.contains("variance_components=refit"));
        assertTrue(log.contains("mixed_fit=laplace-marginal-refit"));
        assertTrue(log.contains("statistic_type=z"));
        assertTrue(log.contains("df_method=asymptotic"));
    }

    @Test
    void pedigreeCliUsesSparseAdditiveRelationshipPrecision() throws Exception {
        Path phenotype = temporaryDirectory.resolve("pedigree-pheno.tsv");
        Path pedigree = temporaryDirectory.resolve("pedigree.tsv");
        Path output = temporaryDirectory.resolve("pedigree-results.tsv");
        StringBuilder observations =
            new StringBuilder("observation\tanimal\ty\tx\n");
        for (int animal = 1; animal <= 11; animal++)
            for (int repeat = 0; repeat < 3; repeat++)
                observations.append('O').append(animal).append('_').append(repeat)
                    .append('\t').append(animal).append('\t')
                    .append(1.0 + animal * 0.18 + repeat * 0.11
                        + ((animal * 7 + repeat * 3) % 5 - 2) * 0.07)
                    .append('\t').append(repeat - 1).append('\n');
        StringBuilder ancestry = new StringBuilder("animal\tsire\tdam\n");
        ancestry.append("1\t0\t0\n2\t0\t0\n");
        for (int animal = 3; animal <= 10; animal++)
            ancestry.append(animal).append('\t')
                .append(Math.max(1, animal - 2)).append('\t')
                .append(Math.max(2, animal - 1)).append('\n');
        Files.writeString(phenotype, observations);
        Files.writeString(pedigree, ancestry);

        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream console = new PrintStream(bytes);
        int status = JLinAlgCli.run(new String[] {
            "--pheno", phenotype.toString(), "--id", "observation",
            "--individual-id", "animal",
            "--formula", "y ~ x + (1|animal)",
            "--pedigree", pedigree.toString(),
            "--pedigree-id", "animal", "--sire-id", "sire",
            "--dam-id", "dam", "--backend", "cpu",
            "--out", output.toString()
        }, console, console);

        assertEquals(0, status);
        assertTrue(bytes.toString().contains(
            "Pedigree members: 11 (file=10, singleton families=1, "
                + "singleton observations=3, "
                + "matching column=animal)"));
        String log = Files.readString(Path.of(output + ".log"));
        assertTrue(log.contains("pedigree_file_members=10"));
        assertTrue(log.contains("pedigree_singletons_added=1"));
        assertTrue(log.contains("pedigree_singleton_observations=3"));
        assertTrue(log.contains("pedigree_members=11"));
        assertTrue(log.contains(
            "pedigree_precision=sparse_additive_relationship_inverse"));
        String manifest = Files.readString(Path.of(output + ".manifest.json"));
        assertTrue(manifest.contains(
            "\"pedigree_file_members\": \"10\""));
        assertTrue(manifest.contains(
            "\"pedigree_singletons_added\": \"1\""));
        assertTrue(manifest.contains(
            "\"pedigree_singleton_observations\": \"3\""));
        assertTrue(manifest.contains(
            "\"pedigree_members\": \"11\""));
        assertTrue(manifest.contains(
            "\"pedigree_precision\": \"sparse_additive_relationship_inverse\""));
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
