/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

/** Result of a balanced one-way ANOVA power calculation. */
public record AnovaPowerResult(
        double groups,
        double observationsPerGroup,
        double betweenVariance,
        double withinVariance,
        double significanceLevel,
        double power) {

    /** Base-R-compatible description of the calculation. */
    public String method() {
        return "Balanced one-way analysis of variance power calculation";
    }

    /** Clarifies that {@link #observationsPerGroup()} is a per-group count. */
    public String note() { return "observationsPerGroup is the number in each group"; }
}
