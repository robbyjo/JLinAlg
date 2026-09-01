/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

/** Nested maximum-likelihood mixed-model likelihood-ratio comparison. */
public record MixedModelComparisonResult(
        double reducedLogLikelihood,
        double fullLogLikelihood,
        double likelihoodRatioStatistic,
        int degreesOfFreedom,
        double pValue) { }
