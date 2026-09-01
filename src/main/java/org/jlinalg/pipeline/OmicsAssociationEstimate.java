/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** One transcript, methylation, protein, or other feature association. */
public record OmicsAssociationEstimate(
        String featureId,
        double beta,
        double standardError,
        double statistic,
        double degreesOfFreedom,
        double pValue,
        double log10PValue,
        double negativeLog10PValue) { }
