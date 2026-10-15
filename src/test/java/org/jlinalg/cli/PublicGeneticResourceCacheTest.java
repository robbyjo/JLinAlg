package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;
import java.security.MessageDigest;
import java.util.HexFormat;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.io.TempDir;

class PublicGeneticResourceCacheTest {
    @TempDir Path directory;
    private static final String SOURCE = "https://example.org/public-reference.tsv";

    @AfterEach void releaseProtectedResource() throws Exception {
        Path file = directory.resolve("reference.tsv");
        if (Files.exists(file)) PublicGeneticResourceCache.releaseVerification(file);
    }

    @Test void concurrentLookupsShareTheVerifiedInstallation() throws Exception {
        Path file = fixture();
        var workers = java.util.concurrent.Executors.newFixedThreadPool(4);
        try {
            var requests = new java.util.ArrayList<java.util.concurrent.Future<PublicGeneticResourceCache.Installed>>();
            for (int i = 0; i < 8; i++) requests.add(workers.submit(() -> installed(file, null)));
            var first = requests.get(0).get();
            for (var request : requests) assertSame(first, request.get());
        } finally {
            workers.shutdownNow();
        }
    }

    @Test void unchangedResourceIsVerifiedOnceAndManifestStillChecked() throws Exception {
        Path file = fixture();
        var first = installed(file, null);
        assertSame(first, installed(file, null));
        Files.writeString(manifest(file), Files.readString(manifest(file))
            .replace(SOURCE, "https://example.org/different.tsv"));
        assertThrows(java.io.IOException.class, () -> installed(file, null));
    }

    @Test void sameSizeChangedContentIsReverified() throws Exception {
        Path file = fixture();
        installed(file, null);
        FileTime before = Files.getLastModifiedTime(file);
        try { Files.writeString(file, "abd"); }
        catch (java.io.IOException protectedFile) {
            assertEquals("abc", Files.readString(file));
            PublicGeneticResourceCache.releaseVerification(file);
            Files.writeString(file, "abd");
        }
        Files.setLastModifiedTime(file, FileTime.fromMillis(before.toMillis() + 5000));
        var failure = assertThrows(java.io.IOException.class,
            () -> installed(file, null));
        assertTrue(failure.getMessage().contains("SHA-256 mismatch"));
    }

    @Test void changedTimestampWithSameContentGetsFreshVerification() throws Exception {
        Path file = fixture();
        var first = installed(file, null);
        FileTime next = FileTime.fromMillis(Files.getLastModifiedTime(file).toMillis() + 5000);
        try { Files.setLastModifiedTime(file, next); }
        catch (java.io.IOException protectedFile) {
            PublicGeneticResourceCache.releaseVerification(file);
            Files.setLastModifiedTime(file, next);
        }
        var second = installed(file, null);
        assertNotSame(first, second);
        assertEquals(first.sha256(), second.sha256());
        assertSame(second, installed(file, null));
    }

    @Test void replacedFileIsNotTrustedEvenWithPreservedSizeAndTimestamp() throws Exception {
        Path file = fixture();
        installed(file, null);
        var attributes = Files.readAttributes(file, BasicFileAttributes.class);
        Path replacement = directory.resolve("replacement.tsv");
        Files.writeString(replacement, "abd");
        // Windows does not expose a BasicFileAttributes file key; use its
        // creation-time identity fallback even on coarse timestamp filesystems.
        if (attributes.fileKey() == null) Files.setAttribute(replacement,
            "basic:creationTime", FileTime.fromMillis(attributes.creationTime().toMillis() + 5000));
        Files.setLastModifiedTime(replacement, attributes.lastModifiedTime());
        try { Files.move(replacement, file, StandardCopyOption.REPLACE_EXISTING); }
        catch (java.io.IOException protectedFile) {
            assertEquals("abc", Files.readString(file));
            PublicGeneticResourceCache.releaseVerification(file);
            Files.move(replacement, file, StandardCopyOption.REPLACE_EXISTING);
        }
        assertThrows(java.io.IOException.class, () -> installed(file, null));
    }

    @Test void changedManifestChecksumAndExpectedChecksumCannotReuseVerification()
            throws Exception {
        Path file = fixture();
        installed(file, null);
        assertThrows(java.io.IOException.class, () -> installed(file,
            new PublicGeneticResourceCache.ExpectedChecksum("MD5", "0".repeat(32))));
        Files.writeString(manifest(file), Files.readString(manifest(file))
            .replace(sha("abc"), "0".repeat(64)));
        assertThrows(java.io.IOException.class, () -> installed(file, null));
    }

    private Path fixture() throws Exception {
        Path file = directory.resolve("reference.tsv");
        Files.writeString(file, "abc");
        StudyPublicEvidence.write(manifest(file), Map.of("schema_version", 1,
            "source", SOURCE, "bytes", 3, "sha256", sha("abc")));
        return file;
    }
    private static Path manifest(Path file) {
        return file.resolveSibling(file.getFileName() + ".jlinalg-resource.yaml");
    }
    private static String sha(String value) throws Exception {
        return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256")
            .digest(value.getBytes(java.nio.charset.StandardCharsets.UTF_8)));
    }
    private static PublicGeneticResourceCache.Installed installed(Path file,
            PublicGeneticResourceCache.ExpectedChecksum expected) throws Exception {
        return PublicGeneticResourceCache.installed(file, manifest(file), SOURCE, expected);
    }
}
