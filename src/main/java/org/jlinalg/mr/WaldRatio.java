/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** A single-variant Wald ratio and first-order standard error. */
public record WaldRatio(
        String variantId,
        double estimate,
        double standardError,
        double statistic,
        double pValue,
        double fStatistic) {
}
