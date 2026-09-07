/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mediation;

/** An estimated mediation path or effect with frequentist uncertainty. */
public record MediationEffect(
        String name,
        double estimate,
        double standardError,
        double statistic,
        double pValue,
        double confidenceLower,
        double confidenceUpper,
        double degreesOfFreedom) {
}
