/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Structured variant-exclusion reasons retained by pipeline audits. */
public enum VariantFilterReason {
    NO_CALLED_SAMPLES,
    TOO_MANY_MISSING,
    BELOW_MINIMUM_MAF,
    ABOVE_MAXIMUM_MAF,
    BELOW_MINIMUM_MAC,
    ABOVE_MAXIMUM_MAC,
    BELOW_IMPUTATION_QUALITY,
    MONOMORPHIC
}
