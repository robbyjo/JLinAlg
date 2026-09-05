/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

/** Immutable confidence interval with its nominal coverage level. */
public record ConfidenceInterval(double lower, double upper, double level) {
    public ConfidenceInterval {
        if (Double.isNaN(lower) || Double.isNaN(upper) || lower > upper) {
            throw new IllegalArgumentException("confidence limits must be ordered numbers");
        }
        if (!(level > 0.0 && level < 1.0) || !Double.isFinite(level)) {
            throw new IllegalArgumentException("confidence level must be in (0, 1)");
        }
    }
}
