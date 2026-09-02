/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

/** Chooses a bounded feature block from current JVM heap headroom. */
final class AdaptiveBlockSizer {
    private static final int MAXIMUM_BLOCK = 8192;
    private static final int MINIMUM_BLOCK = 1;
    private static final long FIXED_RESERVE = 64L * 1024 * 1024;

    private AdaptiveBlockSizer() { }

    static int choose(int samples, int requested) {
        if (requested > 0) return requested;
        Runtime runtime = Runtime.getRuntime();
        long currentlyUsed = runtime.totalMemory() - runtime.freeMemory();
        long headroom = Math.max(0, runtime.maxMemory() - currentlyUsed);
        long budget = Math.max(8L * 1024 * 1024,
            (headroom - Math.min(headroom / 2, FIXED_RESERVE)) / 4);
        long bytesPerFeature = Math.max(128L, 8L * samples * 4L + 256L);
        long calculated = budget / bytesPerFeature;
        return (int) Math.max(MINIMUM_BLOCK,
            Math.min(MAXIMUM_BLOCK, calculated));
    }
}
