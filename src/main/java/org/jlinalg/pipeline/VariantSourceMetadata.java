/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;

/** Metadata available without materializing variant dosage rows. */
public record VariantSourceMetadata(
        Path path, long variantCount, List<String> sampleIds, DataFormat format) {
    public VariantSourceMetadata {
        if (path == null || format == null)
            throw new IllegalArgumentException("path and format are required");
        if (variantCount < -1)
            throw new IllegalArgumentException("variant count must be -1 or nonnegative");
        sampleIds = List.copyOf(sampleIds);
        if (sampleIds.isEmpty())
            throw new IllegalArgumentException("at least one sample is required");
        HashSet<String> seen = new HashSet<>();
        for (String sample : sampleIds) {
            if (sample == null || sample.isBlank() || !seen.add(sample))
                throw new IllegalArgumentException(
                    "sample IDs must be unique and nonblank");
        }
        path = path.toAbsolutePath().normalize();
    }
}
