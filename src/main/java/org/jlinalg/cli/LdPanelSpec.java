/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.nio.file.Path;
import java.util.List;
import org.jlinalg.genetics.LdReferenceLayout;

/** Source mapping and stable identity for one ancestry/cohort LD panel. */
record LdPanelSpec(String id, String ancestry, String sourcePrefix) {
    List<String> sourceFiles() {
        return List.of(sourcePrefix + ".bed", sourcePrefix + ".bim",
            sourcePrefix + ".fam");
    }

    Path targetPrefix(Path referenceDirectory) {
        return LdReferenceLayout.panelPrefix(referenceDirectory, id);
    }

    Path targetFile(Path referenceDirectory, String sourceFile) {
        if (!sourceFiles().contains(sourceFile))
            throw new IllegalArgumentException("file does not belong to panel "
                + id + ": " + sourceFile);
        String suffix = sourceFile.substring(sourceFile.lastIndexOf('.'));
        return Path.of(targetPrefix(referenceDirectory) + suffix);
    }
}
