/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;

/** Re-openable metadata-first source of additive dosage blocks. */
public interface VariantSource {
    VariantSourceMetadata metadata();

    /**
     * Opens a reader whose output samples follow source column indexes in
     * {@code sampleOrder}. Null selects every source sample in source order.
     */
    VariantBlockReader open(int[] sampleOrder) throws IOException;

    default VariantBlockReader open() throws IOException { return open(null); }
}
