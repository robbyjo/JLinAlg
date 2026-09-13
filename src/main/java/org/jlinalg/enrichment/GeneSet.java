/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.enrichment;

import java.util.Set;

/** A named, deduplicated collection in a single gene identifier namespace. */
public record GeneSet(String id, String name, Set<String> genes) {
    public GeneSet {
        if (id == null || id.isBlank() || name == null || genes == null
                || genes.stream().anyMatch(g -> g == null || g.isBlank()))
            throw new IllegalArgumentException("gene set IDs and members must be nonblank");
        genes = Set.copyOf(genes);
    }
}
