/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;

/** One outcome phenotype, optionally grouped into a disease/risk-factor family. */
public record XwasMrOutcome(
        String id,
        String label,
        String category,
        List<SummaryAssociation> associations) {
    /** Copies and validates phenotype metadata and association rows. */
    public XwasMrOutcome {
        if (id == null || id.isBlank())
            throw new IllegalArgumentException("outcome id must be nonblank");
        if (label == null || label.isBlank()) label = id;
        if (category == null) category = "";
        if (associations == null || associations.isEmpty()
                || associations.stream().anyMatch(value -> value == null))
            throw new IllegalArgumentException(
                "outcome associations must be nonempty and contain no nulls");
        associations = List.copyOf(associations);
    }
}
