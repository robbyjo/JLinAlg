/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

/** Reference distribution used for coefficient-level Wald statistics. */
public enum StatisticDistribution {
    /** Student's t distribution, potentially with coefficient-specific DF. */
    STUDENT_T,

    /** Standard normal distribution (a Wald z statistic). */
    STANDARD_NORMAL,
    /** Joint Wald chi-square test. */
    CHI_SQUARE,
    /** Joint finite-denominator Wald F test. */
    F
}
