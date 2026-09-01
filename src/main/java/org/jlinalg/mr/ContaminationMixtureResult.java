/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Grid-profile contamination-mixture causal estimate. */
public record ContaminationMixtureResult(
        MrEstimate estimate, double validInstrumentProbability,
        double logLikelihood, int gridPoints) { }
