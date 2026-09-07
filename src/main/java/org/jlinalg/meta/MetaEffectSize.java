/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** A constructed study-level effect and its first-order sampling error. */
public record MetaEffectSize(String method, double effectSize,
                             double standardError) {
    public MetaEffectSize {
        if (method == null || method.isBlank())
            throw new IllegalArgumentException("effect-size method is required");
        if (!Double.isFinite(effectSize))
            throw new IllegalArgumentException("effect size must be finite");
        if (!(standardError > 0.0) || !Double.isFinite(standardError))
            throw new IllegalArgumentException("standard error must be finite and positive");
    }

    /** Converts this value to a study with the supplied identifier. */
    public MetaStudy study(String name) {
        return new MetaStudy(name, effectSize, standardError);
    }
}
