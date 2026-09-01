/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.ols.OlsResult;

/** Conditional active-set OLS refit after LASSO or elastic-net selection. */
public final class PostSelectionOlsResult {
    private final int[] activePredictorIndices;
    private final boolean includesIntercept;
    private final OlsResult ols;

    PostSelectionOlsResult(
            int[] activePredictorIndices,
            boolean includesIntercept,
            OlsResult ols) {
        this.activePredictorIndices = activePredictorIndices.clone();
        this.includesIntercept = includesIntercept;
        this.ols = ols;
    }

    /** Original zero-based predictor columns included in the refit. */
    public int[] activePredictorIndices() {
        return activePredictorIndices.clone();
    }
    public boolean includesIntercept() { return includesIntercept; }
    public OlsResult ols() { return ols; }
    /**
     * OLS association table in refit-column order. If present, the intercept
     * is first, followed by {@link #activePredictorIndices()} order.
     */
    public AssociationStatistics associationStatistics() {
        return ols.associationStatistics();
    }
}
