/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Objects;

/** One study-level effect estimate and its sampling standard error. */
public record MetaStudy(String name, double effectSize, double standardError) {
    public MetaStudy {
        Objects.requireNonNull(name, "name");
        if (name.isBlank())
            throw new IllegalArgumentException("study name must not be blank");
        if (!Double.isFinite(effectSize))
            throw new IllegalArgumentException("effect size must be finite");
        if (!(standardError > 0.0) || !Double.isFinite(standardError))
            throw new IllegalArgumentException(
                "sampling standard error must be finite and positive");
    }

    /** Sampling variance. */
    public double variance() { return standardError * standardError; }
}
