/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.loess;

/** Controls for one-dimensional direct-surface LOESS. Robust iterations use R's convention. */
public record LoessOptions(
        double span,
        int degree,
        LoessFamily family,
        int robustnessIterations) {
    public LoessOptions {
        if (!(span > 0.0 && span <= 1.0) || !Double.isFinite(span)
                || degree < 0 || degree > 2 || family == null
                || robustnessIterations < 0)
            throw new IllegalArgumentException("invalid LOESS controls");
    }

    /** R-compatible defaults for span, degree, and family. */
    public static LoessOptions defaults() {
        return new LoessOptions(0.75, 2, LoessFamily.GAUSSIAN, 4);
    }

    public LoessOptions withSpan(double value) {
        return new LoessOptions(value, degree, family, robustnessIterations);
    }

    public LoessOptions withDegree(int value) {
        return new LoessOptions(span, value, family, robustnessIterations);
    }

    public LoessOptions withFamily(LoessFamily value) {
        return new LoessOptions(span, degree, value, robustnessIterations);
    }
}
