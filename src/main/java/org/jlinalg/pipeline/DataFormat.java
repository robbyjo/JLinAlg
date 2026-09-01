/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.nio.file.Path;
import java.util.Locale;

/** Supported predictor and variant input formats. */
public enum DataFormat {
    CSV,
    TSV,
    VCF,
    VCF_GZ,
    BCF,
    BGEN;

    /** Infers a format from a conventional filename suffix. */
    public static DataFormat infer(Path path) {
        if (path == null || path.getFileName() == null) {
            throw new IllegalArgumentException("input path is required");
        }
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".vcf.gz") || name.endsWith(".vcf.bgz")) return VCF_GZ;
        if (name.endsWith(".vcf")) return VCF;
        if (name.endsWith(".bcf")) return BCF;
        if (name.endsWith(".bgen")) return BGEN;
        if (name.endsWith(".tsv") || name.endsWith(".txt")) return TSV;
        if (name.endsWith(".csv")) return CSV;
        throw new IllegalArgumentException(
            "cannot infer data format from filename: " + path);
    }
}
