/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

/** Variance estimators for prepared Cox score scans. */
public enum CoxScoreVariance {
    /** Efficient model-based information under the fitted null model. */
    MODEL_BASED,
    /** Cluster-sandwich score variance from caller-supplied cluster IDs. */
    CLUSTER_ROBUST,
    /** Caller-supplied observation score-correlation quadratic form. */
    RELATEDNESS
}
