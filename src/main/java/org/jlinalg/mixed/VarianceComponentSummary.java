/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

/** lme4-style scalar variance-component summary. */
public record VarianceComponentSummary(
        String name, double variance, double standardDeviation) {
    public static VarianceComponentSummary of(String name, double variance) {
        return new VarianceComponentSummary(
            name, variance, Math.sqrt(Math.max(0.0, variance)));
    }
}
