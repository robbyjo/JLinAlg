/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gwas;

/** Missing-genotype handling for association scans. */
public enum GenotypeMissingPolicy {
    /** Reject non-finite genotype dosages. */
    ERROR,
    /** Replace missing dosages by the marker mean before residualization. */
    MEAN_IMPUTE
}
