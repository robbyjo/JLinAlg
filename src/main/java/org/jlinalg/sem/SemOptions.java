/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import org.jlinalg.model.MissingDataPolicy;

/** Objective-evaluation budget, per-case score tolerance and complete-case controls. */
public record SemOptions(int maximumEvaluations, double tolerance,
                         MissingDataPolicy missingDataPolicy) {
    public SemOptions {
        if (maximumEvaluations < 100 || !(tolerance > 0.0) || !Double.isFinite(tolerance)
                || missingDataPolicy == null)
            throw new IllegalArgumentException("invalid SEM options");
    }
    public static SemOptions defaults() {
        return new SemOptions(10_000, 1e-8, MissingDataPolicy.ERROR);
    }
}
