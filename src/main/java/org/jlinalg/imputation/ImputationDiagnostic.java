/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

/** Chain summaries for the originally missing cells of one variable. */
public record ImputationDiagnostic(int variableIndex, VariableType type,
        int observedCount, int missingCount, double[] chainMeans,
        double[] chainVariances, double betweenChainVariance,
        double averageWithinChainVariance) {
    public ImputationDiagnostic {
        chainMeans = chainMeans.clone();
        chainVariances = chainVariances.clone();
    }
    @Override public double[] chainMeans() { return chainMeans.clone(); }
    @Override public double[] chainVariances() { return chainVariances.clone(); }
}
