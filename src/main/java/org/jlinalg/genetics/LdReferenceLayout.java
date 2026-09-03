/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import java.nio.file.Path;
import java.util.Objects;
import java.util.regex.Pattern;

/**
 * Stable on-disk layout shared by every JLinAlg LD reference database.
 *
 * <p>Genotypes use variant-major PLINK 1 BED encoding. A panel prefix points
 * to sibling {@code .bed}, {@code .bim}, and {@code .fam} files.</p>
 */
public final class LdReferenceLayout {
    public static final int FORMAT_VERSION = 1;
    public static final String FORMAT_NAME = "jlinalg-ld-reference";
    public static final String MANIFEST_FILE = "jlinalg-ld-reference.json";
    public static final String GENOTYPE_ENCODING =
        "PLINK_1_BED_VARIANT_MAJOR";
    private static final Pattern PANEL_ID =
        Pattern.compile("[A-Za-z0-9][A-Za-z0-9._-]*");

    private LdReferenceLayout() { }

    /**
     * Returns the canonical PLINK prefix for a panel in a reference directory.
     *
     * @param referenceDirectory installed reference root
     * @param panelId stable ancestry or cohort identifier
     * @return {@code referenceDirectory/panels/panelId/genotypes}
     */
    public static Path panelPrefix(Path referenceDirectory, String panelId) {
        Objects.requireNonNull(referenceDirectory, "referenceDirectory");
        Objects.requireNonNull(panelId, "panelId");
        if (!PANEL_ID.matcher(panelId).matches())
            throw new IllegalArgumentException("invalid LD panel ID: "
                + panelId);
        return referenceDirectory.resolve("panels").resolve(panelId)
            .resolve("genotypes");
    }

    /** Returns the manifest path for an installed reference directory. */
    public static Path manifest(Path referenceDirectory) {
        return Objects.requireNonNull(referenceDirectory,
            "referenceDirectory").resolve(MANIFEST_FILE);
    }
}
