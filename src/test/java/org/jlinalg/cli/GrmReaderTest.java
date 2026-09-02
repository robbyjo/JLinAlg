/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class GrmReaderTest {
    @TempDir Path temporaryDirectory;

    @Test
    void readsDenseMatrixAndOrdersRowsByColumnIds() throws Exception {
        Path path = temporaryDirectory.resolve("relatedness.csv");
        Files.writeString(path,
            "IID,S1,S2,S3\n"
            + "S3,0.2,0.3,1.0\n"
            + "S1,1.0,0.4,0.2\n"
            + "S2,0.4,1.0,0.3\n");

        GrmReader.Loaded loaded = GrmReader.read(path);

        assertEquals("dense", loaded.format());
        assertEquals(java.util.List.of("S1", "S2", "S3"),
            loaded.matrix().sampleIds());
        assertArrayEquals(new double[] {
            1.0, 0.4, 0.2,
            0.4, 1.0, 0.3,
            0.2, 0.3, 1.0
        }, loaded.matrix().relationshipMatrix());
    }

    @Test
    void readsGctaLowerTriangleFromPrefix() throws Exception {
        Path prefix = temporaryDirectory.resolve("cohort");
        Files.writeString(Path.of(prefix + ".grm.id"),
            "F1 S1\nF1 S2\nF2 S3\n");
        float[] triangle = {1.0f, 0.25f, 1.0f, 0.1f, 0.2f, 1.0f};
        ByteBuffer bytes = ByteBuffer.allocate(triangle.length * Float.BYTES)
            .order(ByteOrder.LITTLE_ENDIAN);
        for (float value : triangle) bytes.putFloat(value);
        Files.write(Path.of(prefix + ".grm.bin"), bytes.array());

        GrmReader.Loaded loaded = GrmReader.read(prefix);

        assertEquals("gcta-binary", loaded.format());
        assertEquals(java.util.List.of("S1", "S2", "S3"),
            loaded.matrix().sampleIds());
        assertArrayEquals(new double[] {
            1.0, 0.25, 0.1,
            0.25, 1.0, 0.2,
            0.1, 0.2, 1.0
        }, loaded.matrix().relationshipMatrix(), 1e-7);
    }
}
