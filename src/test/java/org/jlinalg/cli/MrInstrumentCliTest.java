/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.zip.GZIPOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MrInstrumentCliTest {
    @TempDir Path temporaryDirectory;

    @Test
    void searchExpandsBmiAndPrintsDownloadableStudy() {
        AtomicReference<String> query = new AtomicReference<>();
        InstrumentStudy study = study();
        InstrumentCatalog catalog = new StubCatalog(study, new byte[0]) {
            @Override public List<InstrumentStudy> search(String trait, int limit) {
                query.set(trait + ":" + limit);
                return List.of(study);
            }
        };
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream output = new PrintStream(bytes);

        int status = MrInstrumentCli.run(new String[] {
            "search", "--trait", "BMI", "--limit", "7"
        }, output, output, temporaryDirectory, catalog);

        assertEquals(0, status);
        assertEquals("body mass index:7", query.get());
        String text = bytes.toString(StandardCharsets.UTF_8);
        assertTrue(text.contains("GCST000001\tBody mass index"));
        assertTrue(text.contains("expanded trait: BMI -> body mass index"));
    }

    @Test
    void formatMapsUserColumnsFiltersAndRejectsNonSnpRows() throws Exception {
        Path input = temporaryDirectory.resolve("own-results.csv");
        Path output = temporaryDirectory.resolve("mr.tsv");
        Files.writeString(input,
            "variant,estimate,stderr,allele,other,p,frequency,total,gene_symbol\n"
            + "rs1,0.10,0.02,a,g,1e-9,0.2,10000,FTO\n"
            + "rs2,0.20,0.03,c,t,1e-4,0.3,10000,MC4R\n"
            + "rs3,0.30,0.04,AT,a,1e-10,0.1,10000,GENE3\n");
        PrintStream quiet = new PrintStream(OutputStream.nullOutputStream());

        int status = MrInstrumentCli.run(new String[] {
            "format", "--input", input.toString(), "--out", output.toString(),
            "--trait", "BMI", "--p-threshold", "5e-8",
            "--map", "SNP=variant,beta=estimate,se=stderr,"
                + "effect_allele=allele,other_allele=other,pval=p,"
                + "eaf=frequency,samplesize=total,gene=gene_symbol"
        }, quiet, quiet, temporaryDirectory, new StubCatalog(study(), new byte[0]));

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(output);
        assertEquals(2, lines.size());
        assertEquals(String.join("\t", InstrumentTableFormatter.COLUMNS), lines.get(0));
        assertTrue(lines.get(1).startsWith("BMI\trs1\t0.10\t0.02\t0.2\tA\tG\t1e-9"));
        assertTrue(lines.get(1).endsWith("\t10000\tFTO"));
    }

    @Test
    void downloadStreamsGzipAndConvertsOddsRatioToLogEffect() throws Exception {
        String source = "rsid\todds_ratio\tstandard_error\teffect_allele\t"
            + "other_allele\tp_value\teffect_allele_frequency\n"
            + "rs10\t1.25\t0.04\tG\tA\t2e-9\t0.4\n"
            + "rs11\t1.10\t0.03\tC\tT\t2e-3\t0.2\n";
        Path output = temporaryDirectory.resolve("downloaded.tsv");
        ByteArrayOutputStream messages = new ByteArrayOutputStream();
        PrintStream progress = new PrintStream(messages);

        int status = MrInstrumentCli.run(new String[] {
            "download", "--study", "gcst000001", "--out", output.toString()
        }, progress, progress, temporaryDirectory,
            new StubCatalog(study(), gzip(source)));

        assertEquals(0, status);
        List<String> lines = Files.readAllLines(output);
        assertEquals(2, lines.size());
        String[] row = lines.get(1).split("\t", -1);
        assertEquals(Math.log(1.25), Double.parseDouble(row[2]), 1e-15);
        assertEquals("Body mass index", row[0]);
        assertTrue(messages.toString(StandardCharsets.UTF_8)
            .contains("LD clumping"));
    }

    private static InstrumentStudy study() {
        return new InstrumentStudy("GCST000001", "Body mass index",
            List.of("body mass index"), List.of("10000 European"),
            "10,000 European ancestry individuals", 1_000_000,
            URI.create("https://ftp.ebi.ac.uk/example/GCST000001"),
            URI.create("https://creativecommons.org/publicdomain/zero/1.0/"));
    }

    private static byte[] gzip(String text) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (GZIPOutputStream output = new GZIPOutputStream(bytes)) {
            output.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return bytes.toByteArray();
    }

    private static class StubCatalog implements InstrumentCatalog {
        private final InstrumentStudy study;
        private final byte[] content;

        StubCatalog(InstrumentStudy study, byte[] content) {
            this.study = study;
            this.content = content.clone();
        }

        @Override public List<InstrumentStudy> search(String trait, int limit) {
            return List.of(study);
        }

        @Override public InstrumentStudy study(String accession) {
            return study;
        }

        @Override public RemoteSummary openSummaryStatistics(
                InstrumentStudy selected) {
            return new RemoteSummary(
                URI.create("https://ftp.ebi.ac.uk/example/GCST000001.tsv.gz"),
                new ByteArrayInputStream(content));
        }
    }
}
