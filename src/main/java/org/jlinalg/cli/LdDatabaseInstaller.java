/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Duration;
import java.time.OffsetDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import org.jlinalg.genetics.LdReferenceLayout;

/** Streaming, checksum-verifying installer for cataloged LD databases. */
final class LdDatabaseInstaller {
    private static final int TAR_BLOCK_SIZE = 512;
    private static final int BUFFER_SIZE = 1024 * 1024;
    private static final long PROGRESS_INTERVAL = 128L * 1024 * 1024;
    private static final byte[] PLINK_BED_VARIANT_MAJOR =
        {(byte) 0x6c, (byte) 0x1b, (byte) 0x01};

    void install(LdDatabaseSpec database, Path location, PrintStream output)
            throws IOException, InterruptedException {
        if (isInstalled(database, location)) {
            output.println(database.id() + " is already installed in "
                + location);
            return;
        }
        rejectPartialInstallation(database, location);
        prepareDirectories(database, location);

        Path archive = location.resolve("." + database.id() + ".tgz.part");
        boolean complete = false;
        Throwable failure = null;
        try {
            output.println("Downloading " + database.id() + " ("
                + database.downloadSize() + ") to " + location);
            download(database, archive, output);
            verifyChecksum(database, archive);
            output.println("Checksum verified; normalizing reference panels "
                + "to JLinAlg LD format...");
            extract(database, archive, location);
            List<PanelStatistics> statistics = validate(database, location);
            writeManifest(database, location, statistics);
            complete = true;
            output.println("Installed " + database.id() + " in " + location);
            output.println("JLinAlg LD format version "
                + LdReferenceLayout.FORMAT_VERSION + " panel prefixes:");
            for (LdPanelSpec panel : database.panels())
                output.println("  " + panel.id() + ": "
                    + panel.targetPrefix(location));
        } catch (IOException | InterruptedException
                | RuntimeException exception) {
            failure = exception;
            throw exception;
        } finally {
            try {
                cleanup(database, location, archive, complete);
            } catch (IOException cleanupFailure) {
                if (failure != null) failure.addSuppressed(cleanupFailure);
                else throw cleanupFailure;
            }
        }
    }

    private static void prepareDirectories(LdDatabaseSpec database,
            Path location) throws IOException {
        Files.createDirectories(location);
        for (LdPanelSpec panel : database.panels())
            Files.createDirectories(panel.targetPrefix(location).getParent());
    }

    private static void download(LdDatabaseSpec database, Path archive,
            PrintStream output) throws IOException, InterruptedException {
        HttpClient client = HttpClient.newBuilder()
            .connectTimeout(Duration.ofSeconds(30))
            .followRedirects(HttpClient.Redirect.NORMAL)
            .build();
        HttpRequest request = HttpRequest.newBuilder(database.downloadUri())
            .header("User-Agent", "JLinAlg-LD-reference-installer/1")
            .GET()
            .build();
        HttpResponse<InputStream> response = client.send(request,
            HttpResponse.BodyHandlers.ofInputStream());
        if (response.statusCode() < 200 || response.statusCode() >= 300) {
            response.body().close();
            throw new IOException("download server returned HTTP "
                + response.statusCode());
        }
        long expected = response.headers().firstValueAsLong("Content-Length")
            .orElse(-1L);
        try (InputStream input = new BufferedInputStream(response.body());
                OutputStream file = new BufferedOutputStream(
                    Files.newOutputStream(archive,
                        StandardOpenOption.CREATE,
                        StandardOpenOption.TRUNCATE_EXISTING))) {
            byte[] buffer = new byte[BUFFER_SIZE];
            long copied = 0;
            long nextProgress = PROGRESS_INTERVAL;
            for (int count = input.read(buffer); count >= 0;
                    count = input.read(buffer)) {
                if (count == 0) continue;
                file.write(buffer, 0, count);
                copied += count;
                if (copied >= nextProgress) {
                    output.println("  downloaded " + humanBytes(copied)
                        + (expected > 0
                            ? " of " + humanBytes(expected) : ""));
                    nextProgress += PROGRESS_INTERVAL;
                }
            }
        }
    }

    private static void verifyChecksum(LdDatabaseSpec database, Path archive)
            throws IOException {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance(database.checksumAlgorithm());
        } catch (NoSuchAlgorithmException exception) {
            throw new IOException("checksum algorithm is unavailable: "
                + database.checksumAlgorithm(), exception);
        }
        try (InputStream input = new BufferedInputStream(
                Files.newInputStream(archive))) {
            byte[] buffer = new byte[BUFFER_SIZE];
            for (int count = input.read(buffer); count >= 0;
                    count = input.read(buffer)) {
                if (count > 0) digest.update(buffer, 0, count);
            }
        }
        String actual = HexFormat.of().formatHex(digest.digest());
        if (!actual.equalsIgnoreCase(database.checksum()))
            throw new IOException("checksum mismatch: expected "
                + database.checksum() + ", got " + actual);
    }

    static void extract(LdDatabaseSpec database, Path archive,
            Path location) throws IOException {
        Set<String> allowed = Set.copyOf(database.sourceFiles());
        Set<String> extracted = new HashSet<>();
        try (InputStream file = new BufferedInputStream(
                    Files.newInputStream(archive));
                InputStream input = new GZIPInputStream(file, BUFFER_SIZE)) {
            byte[] header = new byte[TAR_BLOCK_SIZE];
            while (readHeader(input, header)) {
                String entryName = tarText(header, 0, 100);
                String prefix = tarText(header, 345, 155);
                if (!prefix.isEmpty()) entryName = prefix + "/" + entryName;
                long size = tarSize(header);
                int type = header[156] & 0xff;
                String fileName = baseName(entryName);
                if ((type == 0 || type == '0') && allowed.contains(fileName)) {
                    if (!extracted.add(fileName))
                        throw new IOException(
                            "archive contains duplicate entry: " + fileName);
                    Path target = database.targetFile(location, fileName);
                    Files.createDirectories(target.getParent());
                    try (OutputStream output = new BufferedOutputStream(
                            Files.newOutputStream(partPath(target),
                                StandardOpenOption.CREATE_NEW))) {
                        copyExactly(input, output, size);
                    }
                } else {
                    copyExactly(input, OutputStream.nullOutputStream(), size);
                }
                skipExactly(input, padding(size));
            }
        }
        if (!extracted.equals(allowed))
            throw new IOException("archive is missing expected entries: "
                + difference(allowed, extracted));
        List<Path> moved = new ArrayList<>();
        try {
            for (Path target : database.installedDataFiles(location)) {
                moveAtomically(partPath(target), target);
                moved.add(target);
            }
        } catch (IOException exception) {
            for (Path target : moved) {
                try {
                    Files.deleteIfExists(target);
                } catch (IOException cleanupFailure) {
                    exception.addSuppressed(cleanupFailure);
                }
            }
            throw exception;
        }
    }

    static List<PanelStatistics> validate(LdDatabaseSpec database,
            Path location) throws IOException {
        List<PanelStatistics> result = new ArrayList<>();
        for (LdPanelSpec panel : database.panels()) {
            Path prefix = panel.targetPrefix(location);
            Path bed = Path.of(prefix + ".bed");
            Path bim = Path.of(prefix + ".bim");
            Path fam = Path.of(prefix + ".fam");
            long samples = countLines(fam);
            long variants = countLines(bim);
            if (samples == 0 || variants == 0)
                throw new IOException("empty PLINK metadata for panel "
                    + panel.id());
            try (InputStream input = Files.newInputStream(bed)) {
                byte[] header = input.readNBytes(3);
                if (!Arrays.equals(header, PLINK_BED_VARIANT_MAJOR))
                    throw new IOException("panel " + panel.id()
                        + " is not variant-major PLINK BED");
            }
            long bytesPerVariant = (samples + 3L) / 4L;
            long expectedSize;
            try {
                expectedSize = Math.addExact(3L,
                    Math.multiplyExact(variants, bytesPerVariant));
            } catch (ArithmeticException exception) {
                throw new IOException("PLINK dimensions overflow for panel "
                    + panel.id(), exception);
            }
            long actualSize = Files.size(bed);
            if (actualSize != expectedSize)
                throw new IOException("PLINK BED size mismatch for panel "
                    + panel.id() + ": expected " + expectedSize
                    + " bytes from " + variants + " variants and " + samples
                    + " samples, found " + actualSize);
            result.add(new PanelStatistics(panel, samples, variants));
        }
        return List.copyOf(result);
    }

    private static long countLines(Path path) throws IOException {
        try (Stream<String> lines = Files.lines(path,
                StandardCharsets.UTF_8)) {
            return lines.filter(line -> !line.isBlank()).count();
        }
    }

    static void writeManifest(LdDatabaseSpec database, Path location,
            List<PanelStatistics> statistics) throws IOException {
        StringBuilder json = new StringBuilder("{\n");
        json.append("  \"format\": \"")
            .append(LdReferenceLayout.FORMAT_NAME).append("\",\n")
            .append("  \"format_version\": ")
            .append(LdReferenceLayout.FORMAT_VERSION).append(",\n")
            .append("  \"created\": \"")
            .append(escape(OffsetDateTime.now().toString())).append("\",\n")
            .append("  \"database_id\": \"")
            .append(escape(database.id())).append("\",\n")
            .append("  \"genome_build\": \"")
            .append(escape(database.genomeBuild())).append("\",\n")
            .append("  \"genotype_encoding\": \"")
            .append(LdReferenceLayout.GENOTYPE_ENCODING).append("\",\n")
            .append("  \"variant_metadata\": \"PLINK_BIM\",\n")
            .append("  \"sample_metadata\": \"PLINK_FAM\",\n")
            .append("  \"source_uri\": \"")
            .append(escape(database.downloadUri().toString()))
            .append("\",\n")
            .append("  \"source_checksum\": {\n")
            .append("    \"algorithm\": \"")
            .append(escape(database.checksumAlgorithm())).append("\",\n")
            .append("    \"value\": \"")
            .append(escape(database.checksum())).append("\"\n")
            .append("  },\n")
            .append("  \"panels\": [\n");
        for (int index = 0; index < statistics.size(); index++) {
            PanelStatistics panel = statistics.get(index);
            String prefix = relativePath(location,
                panel.panel().targetPrefix(location));
            json.append("    {\n")
                .append("      \"id\": \"")
                .append(escape(panel.panel().id())).append("\",\n")
                .append("      \"ancestry\": \"")
                .append(escape(panel.panel().ancestry())).append("\",\n")
                .append("      \"prefix\": \"")
                .append(escape(prefix)).append("\",\n")
                .append("      \"sample_count\": ")
                .append(panel.sampleCount()).append(",\n")
                .append("      \"variant_count\": ")
                .append(panel.variantCount()).append("\n")
                .append("    }");
            if (index + 1 < statistics.size()) json.append(',');
            json.append('\n');
        }
        json.append("  ]\n}\n");
        Path manifest = LdReferenceLayout.manifest(location);
        Path part = partPath(manifest);
        Files.writeString(part, json, StandardCharsets.UTF_8,
            StandardOpenOption.CREATE_NEW);
        moveAtomically(part, manifest);
    }

    private static String relativePath(Path root, Path value) {
        return root.toAbsolutePath().normalize()
            .relativize(value.toAbsolutePath().normalize())
            .toString().replace('\\', '/');
    }

    private static String escape(String value) {
        return value.replace("\\", "\\\\").replace("\"", "\\\"")
            .replace("\r", "\\r").replace("\n", "\\n");
    }

    private static boolean readHeader(InputStream input, byte[] header)
            throws IOException {
        int offset = 0;
        while (offset < header.length) {
            int count = input.read(header, offset, header.length - offset);
            if (count < 0) {
                if (offset == 0) return false;
                throw new EOFException("truncated tar header");
            }
            offset += count;
        }
        for (byte value : header) if (value != 0) return true;
        return false;
    }

    private static String tarText(byte[] header, int offset, int length) {
        int end = offset;
        while (end < offset + length && header[end] != 0) end++;
        return new String(header, offset, end - offset,
            StandardCharsets.UTF_8).strip();
    }

    private static long tarSize(byte[] header) throws IOException {
        String value = tarText(header, 124, 12);
        if (value.isEmpty()) return 0;
        try {
            return Long.parseLong(value, 8);
        } catch (NumberFormatException exception) {
            throw new IOException("invalid tar entry size: " + value,
                exception);
        }
    }

    private static void copyExactly(InputStream input, OutputStream output,
            long bytes) throws IOException {
        byte[] buffer = new byte[BUFFER_SIZE];
        long remaining = bytes;
        while (remaining > 0) {
            int count = input.read(buffer, 0,
                (int) Math.min(buffer.length, remaining));
            if (count < 0) throw new EOFException("truncated tar entry");
            if (count == 0) continue;
            output.write(buffer, 0, count);
            remaining -= count;
        }
    }

    private static void skipExactly(InputStream input, long bytes)
            throws IOException {
        copyExactly(input, OutputStream.nullOutputStream(), bytes);
    }

    private static long padding(long size) {
        return (TAR_BLOCK_SIZE - size % TAR_BLOCK_SIZE) % TAR_BLOCK_SIZE;
    }

    private static String baseName(String name) {
        int separator = Math.max(name.lastIndexOf('/'),
            name.lastIndexOf('\\'));
        return name.substring(separator + 1);
    }

    private static String difference(Set<String> expected,
            Set<String> actual) {
        Set<String> missing = new LinkedHashSet<>(expected);
        missing.removeAll(actual);
        return String.join(", ", missing);
    }

    private static boolean isInstalled(LdDatabaseSpec database,
            Path location) throws IOException {
        Path manifest = LdReferenceLayout.manifest(location);
        if (!Files.isRegularFile(manifest)
                || !database.installedDataFiles(location).stream()
                    .allMatch(Files::isRegularFile))
            return false;
        String json = Files.readString(manifest, StandardCharsets.UTF_8);
        return json.contains("\"format\": \""
                + LdReferenceLayout.FORMAT_NAME + "\"")
            && json.contains("\"format_version\": "
                + LdReferenceLayout.FORMAT_VERSION)
            && json.contains("\"database_id\": \""
                + escape(database.id()) + "\"");
    }

    private static void rejectPartialInstallation(LdDatabaseSpec database,
            Path location) throws IOException {
        List<String> present = installedFiles(database, location).stream()
            .filter(Files::exists)
            .map(path -> relativePath(location, path))
            .toList();
        if (!present.isEmpty())
            throw new IOException("installation directory contains a partial "
                + "or conflicting JLinAlg LD database: "
                + String.join(", ", present));
    }

    private static List<Path> installedFiles(LdDatabaseSpec database,
            Path location) {
        List<Path> result = new ArrayList<>();
        result.add(LdReferenceLayout.manifest(location));
        result.addAll(database.installedDataFiles(location));
        return result;
    }

    private static void cleanup(LdDatabaseSpec database, Path location,
            Path archive, boolean complete) throws IOException {
        List<Path> paths = new ArrayList<>();
        paths.add(archive);
        paths.addAll(installedFiles(database, location).stream()
            .map(LdDatabaseInstaller::partPath).toList());
        if (!complete) paths.addAll(installedFiles(database, location));
        deleteFiles(paths);
    }

    private static void deleteFiles(List<Path> paths) throws IOException {
        IOException failure = null;
        for (Path path : paths) {
            try {
                Files.deleteIfExists(path);
            } catch (IOException exception) {
                if (failure == null) failure = exception;
                else failure.addSuppressed(exception);
            }
        }
        if (failure != null) throw failure;
    }

    private static Path partPath(Path target) {
        return target.resolveSibling("." + target.getFileName() + ".part");
    }

    private static void moveAtomically(Path source, Path destination)
            throws IOException {
        try {
            Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException exception) {
            Files.move(source, destination);
        }
    }

    private static String humanBytes(long bytes) {
        return String.format(java.util.Locale.ROOT, "%.1f MB",
            bytes / (1024.0 * 1024.0));
    }

    record PanelStatistics(LdPanelSpec panel, long sampleCount,
            long variantCount) { }
}
