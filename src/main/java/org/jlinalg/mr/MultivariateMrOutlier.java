/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** One instrument's covariance-aware multivariate MR-PRESSO diagnostic. */
public record MultivariateMrOutlier(
        String variantId, double mahalanobisDistance,
        double empiricalPValue, double bonferroniPValue,
        boolean outlier) { }
