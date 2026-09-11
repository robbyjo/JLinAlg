/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

class AdaptiveBlockSizerTest {
    private static final long GIB = 1024L * 1024 * 1024;

    @Test
    void abundantMemoryQueuesTwoCompleteHighCoreWorkerWaves() {
        int block = AdaptiveBlockSizer.choose(
            5027, 0, 96, 128 * GIB, GIB, GIB / 2);

        assertEquals(49152, block);
        int chunk = AdaptiveBlockSizer.chunkSize(block);
        assertEquals(256, chunk);
        assertEquals(192, block / chunk);
        assertEquals(96,
            AdaptiveBlockSizer.workerCapacity(block, chunk, 96));
    }

    @Test
    void smallerMachinesRetainBoundedBaselineBlock() {
        int block = AdaptiveBlockSizer.choose(
            5027, 0, 8, 128 * GIB, GIB, GIB / 2);

        assertEquals(8192, block);
        int chunk = AdaptiveBlockSizer.chunkSize(block);
        assertEquals(256, chunk);
        assertEquals(8,
            AdaptiveBlockSizer.workerCapacity(block, chunk, 8));
    }

    @Test
    void memoryLimitCapsWorkerCapacity() {
        int block = AdaptiveBlockSizer.choose(
            5027, 0, 96, GIB, GIB, 512L * 1024 * 1024);

        assertEquals(728, block);
        int chunk = AdaptiveBlockSizer.chunkSize(block);
        assertEquals(256, chunk);
        assertEquals(3,
            AdaptiveBlockSizer.workerCapacity(block, chunk, 96));
    }

    @Test
    void intermediateMemoryUsesOneCompleteWorkerWave() {
        int block = AdaptiveBlockSizer.choose(
            5027, 0, 96, 16 * GIB, GIB, GIB);

        assertEquals(24576, block);
        int chunk = AdaptiveBlockSizer.chunkSize(block);
        assertEquals(256, chunk);
        assertEquals(96,
            AdaptiveBlockSizer.workerCapacity(block, chunk, 96));
    }

    @Test
    void baselineRoundsUpToCompleteWorkerWaves() {
        int block = AdaptiveBlockSizer.choose(
            5027, 0, 12, 128 * GIB, GIB, GIB / 2);

        assertEquals(9216, block);
        assertEquals(36, block / AdaptiveBlockSizer.chunkSize(block));
    }

    @Test
    void explicitBlockSizeRemainsAuthoritative() {
        int block = AdaptiveBlockSizer.choose(
            5027, 123, 96, 128 * GIB, GIB, GIB / 2);

        assertEquals(123, block);
        int chunk = AdaptiveBlockSizer.chunkSize(block);
        assertEquals(123, chunk);
        assertEquals(1,
            AdaptiveBlockSizer.workerCapacity(block, chunk, 96));
    }
}
