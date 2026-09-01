/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;

/** Metadata for a row-feature by sample numeric matrix. */
public record NumericMatrixMetadata(
        Path path, long rowCount, List<String> sampleIds, DataFormat format) {
    public NumericMatrixMetadata {
        if (path == null || rowCount < 0 || format == null)
            throw new IllegalArgumentException("matrix metadata is invalid");
        if (format != DataFormat.CSV && format != DataFormat.TSV)
            throw new IllegalArgumentException("numeric matrix must be CSV or TSV");
        sampleIds = List.copyOf(sampleIds);
        HashSet<String> seen = new HashSet<>();
        for (String sample : sampleIds)
            if (sample == null || sample.isBlank() || !seen.add(sample))
                throw new IllegalArgumentException(
                    "sample IDs must be unique and nonblank");
        if (sampleIds.isEmpty())
            throw new IllegalArgumentException("matrix requires samples");
        path = path.toAbsolutePath().normalize();
    }
}
