/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

/** Likelihood used to estimate Gaussian covariance components. */
public enum VarianceEstimation {
    /** Restricted maximum likelihood; preferred for variance estimation. */
    REML,
    /** Profile maximum likelihood; appropriate for fixed-effect model comparison. */
    ML
}
