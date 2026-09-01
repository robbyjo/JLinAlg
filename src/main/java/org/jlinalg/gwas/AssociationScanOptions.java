/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gwas;

import java.util.Objects;

/** Memory, missingness, and parallelism controls for batched marker scans. */
public record AssociationScanOptions(
        int batchSize,
        GenotypeMissingPolicy missingPolicy,
        int parallelism) {
    /** Source-compatible controls using one submitting thread. */
    public AssociationScanOptions(
            int batchSize, GenotypeMissingPolicy missingPolicy) {
        this(batchSize, missingPolicy, 1);
    }

    public AssociationScanOptions {
        if (batchSize < 1) throw new IllegalArgumentException("batchSize must be positive");
        if (parallelism < 1) throw new IllegalArgumentException("parallelism must be positive");
        Objects.requireNonNull(missingPolicy, "missingPolicy");
    }

    public static AssociationScanOptions defaults() {
        return new AssociationScanOptions(
            1024, GenotypeMissingPolicy.MEAN_IMPUTE, 1);
    }
}
