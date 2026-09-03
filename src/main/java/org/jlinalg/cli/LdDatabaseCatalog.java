/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.net.URI;
import java.util.List;

/** Catalog of LD references supported by the command-line installer. */
final class LdDatabaseCatalog {
    private static final List<LdDatabaseSpec> DATABASES = List.of(
        new LdDatabaseSpec(
            "1000g-phase3",
            "1000 Genomes Phase 3, MAF > 0.01; AFR/AMR/EAS/EUR/SAS",
            "GRCh37",
            "approximately 1.6 GB",
            URI.create("https://zenodo.org/records/21422694/files/"
                + "1kg.v3.tgz?download=1"),
            "MD5",
            "75e2350c130f1fbeacf4b1bd840ce0f4",
            List.of(
                new LdPanelSpec("AFR", "African", "AFR"),
                new LdPanelSpec("AMR", "Admixed American", "AMR"),
                new LdPanelSpec("EAS", "East Asian", "EAS"),
                new LdPanelSpec("EUR", "European", "EUR"),
                new LdPanelSpec("SAS", "South Asian", "SAS")))
    );

    private LdDatabaseCatalog() { }

    static List<LdDatabaseSpec> databases() {
        return DATABASES;
    }

    static LdDatabaseSpec find(String id) {
        return DATABASES.stream()
            .filter(database -> database.id().equals(id))
            .findFirst()
            .orElse(null);
    }

    static String choices() {
        StringBuilder result = new StringBuilder("Database choices:\n");
        for (LdDatabaseSpec database : DATABASES) {
            result.append("  ").append(database.id()).append("  ")
                .append(database.description()).append("; ")
                .append(database.genomeBuild()).append("; ")
                .append(database.downloadSize()).append('\n');
        }
        return result.toString().stripTrailing();
    }
}
