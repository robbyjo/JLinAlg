/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Covariance estimator used for coefficient inference. */
public enum GeeCovariance {
    /** Inverse model-based sensitivity matrix. */
    NAIVE,
    /** Cluster sandwich covariance. */
    ROBUST,
    /** Degrees-of-freedom rescaled sandwich covariance. */
    DF_ADJUSTED,
    /** Mancl-DeRouen cluster-leverage corrected sandwich covariance. */
    BIAS_CORRECTED,
    /** Kauermann-Carroll inverse-square-root leverage correction. */
    KAUERMANN_CARROLL,
    /** Fay-Graubard bounded diagonal leverage correction. */
    FAY_GRAUBARD,
    /** Exact delete-one-cluster jackknife covariance. */
    JACKKNIFE
}
