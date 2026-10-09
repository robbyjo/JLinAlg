/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Instant;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import org.jlinalg.study.StudyOperationExecutor.Context;

/** Checksum-verified complete public files shared across study projects. */
final class PublicGeneticResourceCache {
    private static final long GIB = 1024L * 1024 * 1024;

    record ExpectedChecksum(String algorithm, String value) {
        ExpectedChecksum {
            algorithm = algorithm.toUpperCase(Locale.ROOT);
            value = value.toLowerCase(Locale.ROOT);
            if (!algorithm.equals("MD5") && !algorithm.equals("SHA-256"))
                throw new IllegalArgumentException(
                    "Resource checksum must use MD5 or SHA-256");
            int length = algorithm.equals("MD5") ? 32 : 64;
            if (!value.matches("[a-f0-9]{" + length + "}"))
                throw new IllegalArgumentException(
                    "Malformed " + algorithm + " resource checksum");
        }
    }

    record Installed(Path path, String source, long bytes, String sha256,
            Map<String, Object> manifest) { }

    private PublicGeneticResourceCache() { }

    static Installed acquire(Context context, String provider, String release,
            String filename, String source, long maximumBytes,
            ExpectedChecksum expected) throws IOException {
        safe(provider, "provider"); safe(release, "release");
        safe(filename, "filename");
        if (maximumBytes < 1 || maximumBytes > 250 * GIB)
            throw new IllegalArgumentException(
                "Public resource size limit must be in (0,250 GiB]");
        java.net.URI uri = java.net.URI.create(source);
        StudyPublicEvidence.validate(uri);
        Path directory = StudyResourceManager.resourceRoot(context.project(),
            context.study()).resolve("public-genetics").resolve(provider)
            .resolve(release).toAbsolutePath().normalize();
        Path file = directory.resolve(filename).normalize();
        if (!file.getParent().equals(directory))
            throw new IllegalArgumentException("Unsafe public resource path");
        Path manifest = directory.resolve(filename + ".jlinalg-resource.yaml");
        Installed present = installed(file, manifest, source, expected);
        if (present != null) return present;
        Files.createDirectories(directory);
        Path lockPath = directory.resolve(".install.lock");
        try (FileChannel channel = FileChannel.open(lockPath,
                    StandardOpenOption.CREATE, StandardOpenOption.WRITE);
                FileLock lock = channel.lock()) {
            if (!lock.isValid())
                throw new IOException("Could not lock public resource cache: "
                    + directory);
            present = installed(file, manifest, source, expected);
            if (present != null) return present;
            if (Files.exists(file) || Files.exists(manifest))
                throw new IOException("Incomplete public resource installation requires review: " + file);
            String blocked = StudyPublicEvidence.resourceInstallBlocked(
                context.study());
            if (!blocked.isEmpty())
                throw new IOException("Public resource is not installed and download is blocked: "
                    + blocked + " (" + file + ")");
            Path stage = Files.createTempFile(directory, ".download-", ".part");
            Path manifestStage = Files.createTempFile(directory,
                ".manifest-", ".part");
            boolean complete = false;
            try {
                PublicResourceDownloader.Result result;
                try {
                    result = PublicResourceDownloader.downloadComplete(uri,
                        stage, maximumBytes);
                } catch (InterruptedException interrupted) {
                    Thread.currentThread().interrupt();
                    throw new IOException("Public resource installation interrupted",
                        interrupted);
                }
                if (expected != null) {
                    String actual = checksum(stage, expected.algorithm());
                    if (!actual.equals(expected.value()))
                        throw new IOException(expected.algorithm()
                            + " mismatch for public resource " + source);
                }
                String sha = checksum(stage, "SHA-256");
                Map<String, Object> record = new LinkedHashMap<>();
                record.put("schema_version", 1);
                record.put("provider", provider);
                record.put("release", release);
                record.put("source", source);
                record.put("filename", filename);
                record.put("bytes", result.bytes());
                record.put("sha256", sha);
                if (expected != null) {
                    record.put("expected_checksum_algorithm",
                        expected.algorithm());
                    record.put("expected_checksum", expected.value());
                }
                if (result.lastModified() != null)
                    record.put("last_modified", result.lastModified());
                record.put("installed_at", Instant.now().toString());
                record.put("scope", "complete public file; no study variants or regions were submitted");
                StudyPublicEvidence.write(manifestStage, record);
                move(stage, file);
                move(manifestStage, manifest);
                complete = true;
                return new Installed(file, source, result.bytes(), sha,
                    Map.copyOf(record));
            } finally {
                if (!complete) {
                    Files.deleteIfExists(stage);
                    Files.deleteIfExists(manifestStage);
                }
            }
        }
    }

    private static Installed installed(Path file, Path manifestPath,
            String source, ExpectedChecksum expected) throws IOException {
        if (!Files.exists(file) && !Files.exists(manifestPath)) return null;
        if (!Files.isRegularFile(file) || !Files.isRegularFile(manifestPath))
            throw new IOException("Incomplete public resource installation requires review: " + file);
        Map<String, Object> manifest;
        try (var reader = Files.newBufferedReader(manifestPath)) {
            manifest = ProjectConfiguration.mapping(
                ProjectConfiguration.yaml().load(reader),
                manifestPath.toString());
        }
        if (!Integer.valueOf(1).equals(manifest.get("schema_version"))
                || !source.equals(manifest.get("source")))
            throw new IOException("Public resource manifest identity mismatch: "
                + manifestPath);
        long bytes = Files.size(file);
        if (!(manifest.get("bytes") instanceof Number recorded)
                || recorded.longValue() != bytes)
            throw new IOException("Public resource size mismatch: " + file);
        String sha = checksum(file, "SHA-256");
        if (!sha.equals(manifest.get("sha256")))
            throw new IOException("Public resource SHA-256 mismatch: " + file);
        if (expected != null
                && !checksum(file, expected.algorithm()).equals(expected.value()))
            throw new IOException("Public resource " + expected.algorithm()
                + " mismatch: " + file);
        return new Installed(file, source, bytes, sha, Map.copyOf(manifest));
    }

    private static String checksum(Path path, String algorithm)
            throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance(algorithm);
            try (var input = Files.newInputStream(path)) {
                byte[] buffer = new byte[4 * 1024 * 1024];
                for (int count; (count = input.read(buffer)) >= 0;)
                    if (count > 0) digest.update(buffer, 0, count);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException impossible) {
            throw new IllegalStateException(impossible);
        }
    }

    private static void move(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException unsupported) {
            Files.move(source, target);
        }
    }

    private static void safe(String value, String label) {
        if (value == null || !value.matches("[A-Za-z0-9][A-Za-z0-9._-]*"))
            throw new IllegalArgumentException("Unsafe public resource "
                + label + ": " + value);
    }
}
