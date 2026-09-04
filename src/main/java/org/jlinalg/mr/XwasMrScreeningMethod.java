/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Fast primary method used before full xWAS MR diagnostics are computed. */
public enum XwasMrScreeningMethod {
    /** Fixed-effect inverse-variance weighted MR. */
    IVW_FIXED,
    /** Multiplicative-random-effect inverse-variance weighted MR. */
    IVW_MULTIPLICATIVE_RANDOM
}
