/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Immutable row-wise phenotype, transcript, methylation, or protein transform. */
@FunctionalInterface
public interface OmicsTransform {
    /** Returns a transformed copy and preserves missing values as NaN. */
    double[] apply(double[] values);

    default OmicsTransform andThen(OmicsTransform next) {
        if (next == null) throw new IllegalArgumentException("next transform is required");
        return values -> next.apply(apply(values));
    }
}
