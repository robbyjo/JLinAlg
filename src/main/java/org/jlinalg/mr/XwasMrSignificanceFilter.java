/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.Objects;
import org.jlinalg.inference.PValueScale;

/** An explicit, scale-aware significance rule for an xWAS MR screen. */
public record XwasMrSignificanceFilter(PValueScale scale, double threshold) {
    /** Validates the threshold in the direction implied by its scale. */
    public XwasMrSignificanceFilter {
        Objects.requireNonNull(scale, "scale");
        if (!Double.isFinite(threshold))
            throw new IllegalArgumentException("threshold must be finite");
        switch (scale) {
            case REGULAR -> {
                if (threshold < 0.0 || threshold > 1.0)
                    throw new IllegalArgumentException(
                        "regular p-value threshold must lie in [0, 1]");
            }
            case LOG10 -> {
                if (threshold > 0.0)
                    throw new IllegalArgumentException(
                        "log10(p) threshold must be nonpositive");
            }
            case NEGATIVE_LOG10 -> {
                if (threshold < 0.0)
                    throw new IllegalArgumentException(
                        "-log10(p) threshold must be nonnegative");
            }
        }
    }

    /** Selects results with {@code p <= maximum}. */
    public static XwasMrSignificanceFilter pValueAtMost(double maximum) {
        return new XwasMrSignificanceFilter(PValueScale.REGULAR, maximum);
    }

    /** Selects results with {@code log10(p) <= maximum}. */
    public static XwasMrSignificanceFilter log10PAtMost(double maximum) {
        return new XwasMrSignificanceFilter(PValueScale.LOG10, maximum);
    }

    /** Selects results with {@code -log10(p) >= minimum}. */
    public static XwasMrSignificanceFilter negativeLog10PAtLeast(
            double minimum) {
        return new XwasMrSignificanceFilter(
            PValueScale.NEGATIVE_LOG10, minimum);
    }

    /** Returns the p-value represented on this filter's scale. */
    public double value(double pValue) {
        validatePValue(pValue);
        return switch (scale) {
            case REGULAR -> pValue;
            case LOG10 -> pValue == 0.0
                ? Double.NEGATIVE_INFINITY : Math.log10(pValue);
            case NEGATIVE_LOG10 -> pValue == 0.0
                ? Double.POSITIVE_INFINITY : -Math.log10(pValue);
        };
    }

    /** Tests the scale-specific inclusive significance boundary. */
    public boolean includes(double pValue) {
        double transformed = value(pValue);
        return scale == PValueScale.NEGATIVE_LOG10
            ? transformed >= threshold : transformed <= threshold;
    }

    private static void validatePValue(double pValue) {
        if (!Double.isFinite(pValue) || pValue < 0.0 || pValue > 1.0)
            throw new IllegalArgumentException("p-value must lie in [0, 1]");
    }
}
