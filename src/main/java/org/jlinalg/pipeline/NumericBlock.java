/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.List;

/** Ordered bounded block of numeric feature rows. */
public record NumericBlock(long firstSourceIndex, List<NumericRow> rows) {
    public NumericBlock {
        if (firstSourceIndex < 0)
            throw new IllegalArgumentException("first source index cannot be negative");
        rows = List.copyOf(rows);
        if (rows.isEmpty()) throw new IllegalArgumentException("block cannot be empty");
    }
}
