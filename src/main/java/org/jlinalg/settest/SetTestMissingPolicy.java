/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/** Missing-dosage handling performed after aligned-cohort frequency QC. */
public enum SetTestMissingPolicy {
    ERROR,
    MEAN_IMPUTE,
    ZERO
}
