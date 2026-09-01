/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

import java.util.Objects;

/** One pair retained from a genomic-relationship threshold query. */
public record RelatednessPair(
        String firstId,
        String secondId,
        double relationship,
        double kinshipCoefficient) {

    public RelatednessPair {
        Objects.requireNonNull(firstId, "firstId");
        Objects.requireNonNull(secondId, "secondId");
        if (!Double.isFinite(relationship)
                || !Double.isFinite(kinshipCoefficient))
            throw new IllegalArgumentException(
                "relatedness estimates must be finite");
    }
}
