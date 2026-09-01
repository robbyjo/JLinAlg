/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.ols;

/** Behavior when an OLS design matrix is not full column rank. */
public enum RankDeficiencyStrategy {
    /** Reject the fit because individual coefficients are not identified. */
    ERROR,
    /** Return the Moore-Penrose minimum-norm solution using SVD. */
    MINIMUM_NORM
}
