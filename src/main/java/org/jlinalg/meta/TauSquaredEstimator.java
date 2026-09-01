/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** Estimator for between-study heterogeneity variance. */
public enum TauSquaredEstimator {
    /** Restricted maximum likelihood; the default for random-effects fits. */
    REML,
    /** DerSimonian-Laird generalized method of moments. */
    DERSIMONIAN_LAIRD,
    /** Paule-Mandel generalized Q estimator. */
    PAULE_MANDEL
}
