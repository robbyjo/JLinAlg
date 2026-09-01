/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Structured reasons why an exposure variant was not harmonized. */
public enum HarmonizationExclusionReason {
    DUPLICATE_EXPOSURE,
    DUPLICATE_OUTCOME,
    MISSING_OUTCOME,
    ALLELE_MISMATCH,
    PALINDROMIC_AMBIGUOUS,
    FREQUENCY_MISMATCH,
    ZERO_EXPOSURE_EFFECT
}
