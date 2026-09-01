/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

/** Factory for format-neutral variant sources. */
public final class VariantSources {
    private VariantSources() { }

    public static VariantSource open(Path path) throws IOException {
        return open(path, DataFormat.infer(path), null);
    }

    public static VariantSource open(Path path, DataFormat format)
            throws IOException {
        return open(path, format, null);
    }

    /**
     * Opens a source. External sample IDs are used only for BGEN files that do
     * not contain an embedded sample identifier block.
     */
    public static VariantSource open(
            Path path, DataFormat format, List<String> externalSampleIds)
            throws IOException {
        return switch (format) {
            case CSV, TSV -> new DelimitedVariantSource(path, format);
            case VCF, VCF_GZ, BCF -> new VcfVariantSource(path, format);
            case BGEN -> new BgenVariantSource(path, externalSampleIds);
        };
    }
}
