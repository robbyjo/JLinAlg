/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

/** One input-ordered coefficient estimate from an association batch. */
public record AssociationEstimate(
        String name,
        double beta,
        double standardError,
        double statistic,
        double pValue,
        double log10PValue,
        double negativeLog10PValue,
        double degreesOfFreedom) { }
