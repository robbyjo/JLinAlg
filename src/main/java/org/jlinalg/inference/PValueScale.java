/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.inference;

/** Numeric representation requested for a p-value result vector. */
public enum PValueScale {
    /** Two-sided p in the interval [0, 1]. */
    REGULAR,
    /** Base-10 logarithm of p, normally non-positive. */
    LOG10,
    /** Negative base-10 logarithm of p, conventional for GWAS plots. */
    NEGATIVE_LOG10
}
