/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

/** Covariance estimator for individual-level two-stage least squares. */
public enum InstrumentalVariableCovariance {
    /** Homoskedastic IV covariance, with residual Student-t inference. */
    HOMOSKEDASTIC,

    /** Eicker-White HC0 covariance, with asymptotic normal inference. */
    HC0,

    /** HC0 multiplied by {@code n / (n - p)}, with asymptotic normal inference. */
    HC1,

    /** Cluster score sandwich without a finite-sample covariance multiplier. */
    CLUSTER_CR0,

    /** Cluster score sandwich with the conventional CR1 finite-sample multiplier. */
    CLUSTER_CR1;

    /** Whether this estimator requires one cluster identifier per observation. */
    public boolean requiresClusters() {
        return this == CLUSTER_CR0 || this == CLUSTER_CR1;
    }

    /** Whether inference is robust to observation-level heteroskedasticity. */
    public boolean isRobust() {
        return this != HOMOSKEDASTIC;
    }
}
