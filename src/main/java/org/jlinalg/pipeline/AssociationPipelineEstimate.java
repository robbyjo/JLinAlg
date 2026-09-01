/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Variant metadata, QC, and one coefficient-level association result. */
public record AssociationPipelineEstimate(
        VariantRecord variant,
        VariantStatistics variantStatistics,
        double beta,
        double standardError,
        double statistic,
        double degreesOfFreedom,
        double pValue,
        double log10PValue,
        double negativeLog10PValue) { }
