/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

import java.util.ArrayList;
import java.util.List;

/** Completed datasets and chain diagnostics from reproducible MICE streams. */
public final class MiceResult {
    private final List<double[][]> datasets;
    private final List<ImputationDiagnostic> diagnostics;
    private final MiceOptions options;

    MiceResult(List<double[][]> datasets,
            List<ImputationDiagnostic> diagnostics, MiceOptions options) {
        this.datasets = new ArrayList<>();
        for (double[][] dataset : datasets) this.datasets.add(copy(dataset));
        this.diagnostics = List.copyOf(diagnostics);
        this.options = options;
    }

    public List<double[][]> datasets() {
        List<double[][]> result = new ArrayList<>();
        for (double[][] dataset : datasets) result.add(copy(dataset));
        return List.copyOf(result);
    }
    public List<ImputationDiagnostic> diagnostics() { return diagnostics; }
    public MiceOptions options() { return options; }

    private static double[][] copy(double[][] source) {
        double[][] result = new double[source.length][];
        for (int row = 0; row < source.length; row++) result[row] = source[row].clone();
        return result;
    }
}
