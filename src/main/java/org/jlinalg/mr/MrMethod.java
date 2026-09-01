/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Summary-data Mendelian randomization estimators. */
public enum MrMethod {
    WALD_RATIO,
    IVW_FIXED,
    IVW_MULTIPLICATIVE_RANDOM,
    IVW_GENERALIZED_FIXED,
    IVW_GENERALIZED_MULTIPLICATIVE_RANDOM,
    MR_EGGER,
    MR_EGGER_GENERALIZED,
    WEIGHTED_MEDIAN,
    MR_RAPS,
    CONTAMINATION_MIXTURE,
    MULTIVARIABLE_IVW,
    MULTIVARIABLE_EGGER,
    OVERLAP_AWARE_IVW
}
