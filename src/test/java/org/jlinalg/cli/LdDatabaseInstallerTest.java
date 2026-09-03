/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.OutputStream;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;
import org.jlinalg.genetics.LdReferenceLayout;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class LdDatabaseInstallerTest {
    @TempDir Path temporaryDirectory;

    @Test
    void extractsOnlyCatalogedFilesFromTarGzip() throws Exception {
        Path archive = temporaryDirectory.resolve("reference.tgz");
        Path destination = temporaryDirectory.resolve("installed");
        Map<String, byte[]> entries = new LinkedHashMap<>();
        entries.put("bundle/EUR.bed",
            new byte[] {(byte) 0x6c, (byte) 0x1b, 0x01, 0x00});
        entries.put("bundle/EUR.bim",
            "1\trs1\t0\t100\tA\tG\n".getBytes(StandardCharsets.UTF_8));
        entries.put("bundle/EUR.fam",
            "F1 I1 0 0 0 -9\n".getBytes(StandardCharsets.UTF_8));
        entries.put("bundle/README.txt",
            "ignored".getBytes(StandardCharsets.UTF_8));
        writeTarGzip(archive, entries);
        Files.createDirectories(destination);
        LdDatabaseSpec database = new LdDatabaseSpec(
            "test", "test", "GRCh37", "small",
            URI.create("https://example.invalid/reference.tgz"),
            "MD5", "unused",
            List.of(new LdPanelSpec("EUR", "European", "EUR")));

        LdDatabaseInstaller.extract(database, archive, destination);
        List<LdDatabaseInstaller.PanelStatistics> statistics =
            LdDatabaseInstaller.validate(database, destination);
        LdDatabaseInstaller.writeManifest(
            database, destination, statistics);

        Path prefix = LdReferenceLayout.panelPrefix(destination, "EUR");
        assertEquals(4, Files.size(Path.of(prefix + ".bed")));
        assertEquals("1\trs1\t0\t100\tA\tG\n",
            Files.readString(Path.of(prefix + ".bim")));
        assertEquals(1, statistics.get(0).sampleCount());
        assertEquals(1, statistics.get(0).variantCount());
        String manifest = Files.readString(
            LdReferenceLayout.manifest(destination));
        assertTrue(manifest.contains(
            "\"genotype_encoding\": \"PLINK_1_BED_VARIANT_MAJOR\""));
        assertTrue(manifest.contains(
            "\"prefix\": \"panels/EUR/genotypes\""));
        assertFalse(Files.exists(destination.resolve("README.txt")));
        assertFalse(Files.exists(
            prefix.getParent().resolve(".genotypes.bed.part")));
    }

    private static void writeTarGzip(Path path, Map<String, byte[]> entries)
            throws IOException {
        try (OutputStream file = Files.newOutputStream(path);
                OutputStream output = new GZIPOutputStream(file)) {
            for (Map.Entry<String, byte[]> entry : entries.entrySet()) {
                byte[] header = new byte[512];
                put(header, 0, 100, entry.getKey());
                put(header, 100, 8, "0000644");
                put(header, 108, 8, "0000000");
                put(header, 116, 8, "0000000");
                put(header, 124, 12,
                    String.format("%011o", entry.getValue().length));
                put(header, 136, 12, "00000000000");
                header[156] = '0';
                put(header, 257, 6, "ustar");
                output.write(header);
                output.write(entry.getValue());
                int padding = (512 - entry.getValue().length % 512) % 512;
                output.write(new byte[padding]);
            }
            output.write(new byte[1024]);
        }
    }

    private static void put(byte[] target, int offset, int length,
            String value) {
        byte[] source = value.getBytes(StandardCharsets.US_ASCII);
        System.arraycopy(source, 0, target, offset,
            Math.min(source.length, length));
    }
}
