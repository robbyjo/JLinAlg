/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

/** Chooses a bounded feature block from JVM heap headroom and scan threads. */
final class AdaptiveBlockSizer {
    private static final int BASELINE_BLOCK_SIZE = 8192;
    private static final int FEATURES_PER_CHUNK = 256;
    private static final int MINIMUM_BLOCK = 1;
    private static final long FIXED_RESERVE = 64L * 1024 * 1024;

    private AdaptiveBlockSizer() { }

    static int choose(int samples, int requested, int threads) {
        Runtime runtime = Runtime.getRuntime();
        return choose(samples, requested, threads, runtime.maxMemory(),
            runtime.totalMemory(), runtime.freeMemory());
    }

    static int choose(
            int samples, int requested, int threads,
            long maximumMemory, long allocatedMemory, long freeMemory) {
        if (requested > 0) return requested;
        if (samples < 1 || threads < 1)
            throw new IllegalArgumentException(
                "samples and threads must be positive");
        long currentlyUsed = Math.max(0, allocatedMemory - freeMemory);
        long headroom = Math.max(0, maximumMemory - currentlyUsed);
        long budget = Math.max(8L * 1024 * 1024,
            (headroom - Math.min(headroom / 2, FIXED_RESERVE)) / 4);
        long bytesPerFeature = Math.max(128L, 8L * samples * 4L + 256L);
        long calculated = budget / bytesPerFeature;
        long threadTarget = Math.max(BASELINE_BLOCK_SIZE,
            Math.min(Integer.MAX_VALUE,
                (long) threads * FEATURES_PER_CHUNK));
        return (int) Math.max(MINIMUM_BLOCK,
            Math.min(threadTarget, calculated));
    }

    static int chunkSize(int blockSize) {
        if (blockSize < 1)
            throw new IllegalArgumentException(
                "block size must be positive");
        return Math.min(FEATURES_PER_CHUNK, blockSize);
    }

    static int workerCapacity(int blockSize, int chunkSize, int threads) {
        if (blockSize < 1 || chunkSize < 1 || threads < 1)
            throw new IllegalArgumentException(
                "block size, chunk size, and threads must be positive");
        long chunks = ((long) blockSize + chunkSize - 1) / chunkSize;
        return (int) Math.min(threads, chunks);
    }
}
