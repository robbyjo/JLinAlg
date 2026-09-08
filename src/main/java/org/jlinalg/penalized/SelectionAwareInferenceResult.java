/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import org.jlinalg.ols.OlsResult;

/** Deterministic data-splitting inference after penalized selection. */
public record SelectionAwareInferenceResult(PenalizedRegressionResult selectionFit,
        int[] activePredictorIndices, int selectionObservations,
        int inferenceObservations, OlsResult inferenceFit) {
    public SelectionAwareInferenceResult { activePredictorIndices = activePredictorIndices.clone(); }
    public int[] activePredictorIndices() { return activePredictorIndices.clone(); }
}
