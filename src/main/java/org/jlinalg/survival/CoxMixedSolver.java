/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Numerical equation strategy used for a Gaussian-frailty Cox fit. */
public enum CoxMixedSolver {
    /** Full dense random-effect information and factorization. */
    DENSE,
    /** Sparse precision with a diagonal random-information approximation. */
    SPARSE_PRECISION
}
