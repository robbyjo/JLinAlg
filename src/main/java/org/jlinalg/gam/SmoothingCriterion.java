/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

/** Criteria available for direct Gaussian smoothing-parameter selection. */
public enum SmoothingCriterion {
    /** Generalized cross-validation with unknown residual scale. */
    GCV,
    /** Unbiased risk estimation with a caller-supplied residual scale. */
    UBRE,
    /** Akaike information criterion using effective degrees of freedom. */
    AIC
}
