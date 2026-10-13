/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

class PublicResourceDownloaderTest {
    @Test
    void boundsCompleteFileWorkerCountByChunksAndConcurrencyLimit() {
        assertEquals(1, PublicResourceDownloader.completeDownloadWorkers(1));
        assertEquals(1, PublicResourceDownloader.completeDownloadWorkers(
            32L * 1024 * 1024));
        assertEquals(2, PublicResourceDownloader.completeDownloadWorkers(
            32L * 1024 * 1024 + 1));
        assertEquals(32, PublicResourceDownloader.completeDownloadWorkers(
            100L * 1024 * 1024 * 1024));
        assertThrows(IllegalArgumentException.class,
            () -> PublicResourceDownloader.completeDownloadWorkers(0));
    }
}
