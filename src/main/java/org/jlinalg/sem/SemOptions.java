/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import org.jlinalg.model.MissingDataPolicy;

/** Optimization and complete-case controls for SEM. */
public record SemOptions(int maximumEvaluations, double tolerance,
                         MissingDataPolicy missingDataPolicy) {
    public SemOptions {
        if (maximumEvaluations < 100 || !(tolerance > 0.0)
                || missingDataPolicy == null)
            throw new IllegalArgumentException("invalid SEM options");
    }
    public static SemOptions defaults() {
        return new SemOptions(10_000, 1e-8, MissingDataPolicy.ERROR);
    }
}
