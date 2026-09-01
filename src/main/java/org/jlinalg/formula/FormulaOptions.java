/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.Objects;

/** Fixed-effect formula compilation controls. */
public record FormulaOptions(
        ContrastCoding contrastCoding,
        String weightColumn) {
    public FormulaOptions {
        Objects.requireNonNull(contrastCoding, "contrastCoding");
        if (weightColumn != null && weightColumn.isBlank()) {
            throw new IllegalArgumentException("weightColumn must not be blank");
        }
    }

    public static FormulaOptions defaults() {
        return new FormulaOptions(ContrastCoding.TREATMENT, null);
    }
}
