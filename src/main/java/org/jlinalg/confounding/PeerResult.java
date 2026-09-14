/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.List;

/** Variational PEER posterior means and convergence diagnostics. */
public record PeerResult(
        double[][] factors,
        double[][] loadings,
        double[][] residuals,
        double[] factorPrecisions,
        double[] featureNoisePrecisions,
        int iterations,
        boolean converged,
        List<Double> evidenceBounds,
        List<Double> residualVariances) {
    public PeerResult {
        evidenceBounds = List.copyOf(evidenceBounds);
        residualVariances = List.copyOf(residualVariances);
    }
}
