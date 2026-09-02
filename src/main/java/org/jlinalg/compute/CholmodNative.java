/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.compute;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/** Loads the optional SuiteSparse/CHOLMOD JNI bridge. */
final class CholmodNative {
    private static final LoadResult LOAD = load();

    private CholmodNative() { }

    static boolean available() { return LOAD.available(); }
    static String description() { return LOAD.description(); }

    static void requireAvailable() {
        if (!LOAD.available()) throw new IllegalStateException(
            "CHOLMOD is unavailable: " + LOAD.description());
    }

    static native long create(
        int dimension, int[] rowStarts, int[] columnIndices, double[] values,
        boolean lower, boolean naturalOrdering);
    static native void refactor(long handle, double[] values);
    static native double logDeterminant(long handle);
    static native int factorNonzeroCount(long handle);
    static native int[] permutation(long handle);
    static native void solveInPlace(long handle, double[] rightHandSide,
        int rightHandSides);
    static native void destroy(long handle);
    private static native String version();

    private static LoadResult load() {
        List<String> failures = new ArrayList<>();
        String explicit = System.getProperty("jlinalg.cholmod.library");
        if (explicit != null && !explicit.isBlank()) {
            LoadResult result = loadPath(Path.of(explicit), failures);
            if (result != null) return result;
        }
        try {
            System.loadLibrary("jlinalg_cholmod");
            return new LoadResult(true, version());
        } catch (LinkageError failure) {
            failures.add(failure.getMessage());
        }
        for (Path candidate : developmentCandidates()) {
            if (!Files.isRegularFile(candidate)) continue;
            LoadResult result = loadPath(candidate, failures);
            if (result != null) return result;
        }
        String resource = resourcePath();
        if (resource != null) {
            try (InputStream input = CholmodNative.class.getResourceAsStream(resource)) {
                if (input != null) {
                    String suffix = System.mapLibraryName("jlinalg_cholmod");
                    Path extracted = Files.createTempFile("jlinalg-cholmod-", suffix);
                    Files.copy(input, extracted, StandardCopyOption.REPLACE_EXISTING);
                    extracted.toFile().deleteOnExit();
                    LoadResult result = loadPath(extracted, failures);
                    if (result != null) return result;
                }
            } catch (IOException | LinkageError failure) {
                failures.add(failure.getMessage());
            }
        }
        String message = failures.stream().filter(value -> value != null)
            .reduce((left, right) -> left + "; " + right)
            .orElse("native library was not found");
        return new LoadResult(false, message);
    }

    private static LoadResult loadPath(Path path, List<String> failures) {
        try {
            System.load(path.toAbsolutePath().normalize().toString());
            return new LoadResult(true, version());
        } catch (LinkageError failure) {
            failures.add(path + ": " + failure.getMessage());
            return null;
        }
    }

    private static List<Path> developmentCandidates() {
        String mapped = System.mapLibraryName("jlinalg_cholmod");
        return List.of(
            Path.of("build", "native", "cholmod", "bin", "Release", mapped),
            Path.of("build", "native", "cholmod", "bin", mapped),
            Path.of("build", "native", "cholmod", "Release", mapped));
    }

    private static String resourcePath() {
        String os = System.getProperty("os.name", "").toLowerCase(Locale.ROOT);
        String arch = System.getProperty("os.arch", "").toLowerCase(Locale.ROOT);
        String platform;
        if (os.contains("win")) platform = "windows";
        else if (os.contains("mac")) platform = "macos";
        else if (os.contains("linux")) platform = "linux";
        else return null;
        String machine = arch.equals("amd64") || arch.equals("x86_64")
            ? "x86_64" : arch.equals("aarch64") || arch.equals("arm64")
                ? "arm64" : arch;
        return "/native/" + platform + "-" + machine + "/"
            + System.mapLibraryName("jlinalg_cholmod");
    }

    private record LoadResult(boolean available, String description) { }
}
