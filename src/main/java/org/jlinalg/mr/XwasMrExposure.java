/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;

/** One xWAS exposure and its already-clumped variant associations. */
public record XwasMrExposure(
        String id,
        String label,
        List<SummaryAssociation> clumpedInstruments) {
    /** Copies and validates exposure metadata and instrument rows. */
    public XwasMrExposure {
        if (id == null || id.isBlank())
            throw new IllegalArgumentException("exposure id must be nonblank");
        if (label == null || label.isBlank()) label = id;
        if (clumpedInstruments == null || clumpedInstruments.isEmpty()
                || clumpedInstruments.stream().anyMatch(value -> value == null))
            throw new IllegalArgumentException(
                "clumped instruments must be nonempty and contain no nulls");
        clumpedInstruments = List.copyOf(clumpedInstruments);
    }
}
