/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** One transcript, methylation site, protein, metabolite, or numeric feature. */
public record NumericRow(String id, double[] values) {
    public NumericRow {
        if (id == null || id.isBlank() || values == null || values.length == 0)
            throw new IllegalArgumentException("row ID and values are required");
        id = id.trim();
        values = values.clone();
    }
    @Override public double[] values() { return values.clone(); }
    double[] valuesView() { return values; }
}
