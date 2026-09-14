/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.confounding;

import java.util.List;

/** Standard or automatically selected surrogate-variable result. */
public record SvaResult(
        double[][] surrogateVariables,
        double[] probabilityHeterogeneity,
        double[] probabilityBiological,
        double[] featureWeights,
        double[] varianceExplained,
        double[][] adjusted,
        int factorCount,
        int iterations,
        boolean converged,
        double finalWeightMrse,
        List<FactorSelection> selectionPath) {
    public SvaResult {
        selectionPath = List.copyOf(selectionPath);
    }

    /** One evaluated factor count in AutoSVA's preserved-signal search. */
    public record FactorSelection(int factors, double fRatio, boolean succeeded) { }
}
