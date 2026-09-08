/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** Cluster sandwich corrections following clubSandwich naming. */
public enum MetaRobustCorrection {
    CR0,
    /** Multiply CR0 by G/(G-1). */
    CR1,
    /** Multiply CR0 by G/(G-p), as metafor robust(adjust=TRUE). */
    CR1P,
    /** Multiply CR0 by G(N-1)/((G-1)(N-p)). */
    CR1S,
    /** Bias-reduced linearization under the fitted covariance target, with Satterthwaite t DF. */
    CR2
}
