/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.net.URI;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Immutable metadata for a downloadable linkage-disequilibrium reference. */
record LdDatabaseSpec(
        String id,
        String description,
        String genomeBuild,
        String downloadSize,
        URI downloadUri,
        String checksumAlgorithm,
        String checksum,
        List<LdPanelSpec> panels) {
    LdDatabaseSpec {
        panels = List.copyOf(panels);
    }

    List<String> sourceFiles() {
        return panels.stream().flatMap(panel -> panel.sourceFiles().stream())
            .toList();
    }

    Path targetFile(Path referenceDirectory, String sourceFile) {
        for (LdPanelSpec panel : panels)
            if (panel.sourceFiles().contains(sourceFile))
                return panel.targetFile(referenceDirectory, sourceFile);
        throw new IllegalArgumentException(
            "unexpected LD reference source file: " + sourceFile);
    }

    List<Path> installedDataFiles(Path referenceDirectory) {
        List<Path> result = new ArrayList<>();
        for (LdPanelSpec panel : panels)
            for (String sourceFile : panel.sourceFiles())
                result.add(panel.targetFile(referenceDirectory, sourceFile));
        return result;
    }
}
