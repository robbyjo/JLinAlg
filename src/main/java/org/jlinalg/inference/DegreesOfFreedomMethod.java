/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

/** Denominator degrees-of-freedom policy for fixed-effect Wald tests. */
public enum DegreesOfFreedomMethod {
    /** Conventional residual DF, normally {@code observations - rank(X)}. */
    RESIDUAL,

    /** Residual DF based on the effective (hat-matrix trace) model dimension. */
    EFFECTIVE_RESIDUAL,

    /**
     * Fast mixed-model approximation {@code observations - rank(X) - 1}.
     * This is the default for REML and PQL working-REML fits.
     */
    RESIDUAL_APPROXIMATION,

    /**
     * Coefficient-specific Satterthwaite degrees of freedom obtained by a
     * first-order delta method over the estimated variance components.
     */
    SATTERTHWAITE,

    /**
     * Kenward-Roger small-sample covariance adjustment and coefficient-specific
     * denominator degrees of freedom for Gaussian REML fits.
     */
    KENWARD_ROGER,

    /** Asymptotic normal inference without finite denominator degrees of freedom. */
    ASYMPTOTIC
}
