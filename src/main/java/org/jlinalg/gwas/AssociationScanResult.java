/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gwas;

import java.util.List;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.PValueScale;
import org.jlinalg.reml.RemlResult;

/** Marker-ordered effect estimates from a variance-component null model. */
public final class AssociationScanResult {
    private final List<String> markerNames;
    private final AssociationStatistics statistics;
    private final RemlResult nullModel;

    AssociationScanResult(
            List<String> markerNames,
            AssociationStatistics statistics,
            RemlResult nullModel) {
        this.markerNames = List.copyOf(markerNames);
        this.statistics = statistics;
        this.nullModel = nullModel;
    }

    public List<String> markerNames() { return markerNames; }
    public AssociationStatistics statistics() { return statistics; }
    public double[] beta() { return statistics.beta(); }
    public double[] effectSizes() { return statistics.effectSizes(); }
    public double[] standardErrors() { return statistics.standardErrors(); }
    public double[] tStatistics() { return statistics.statistics(); }
    public double[] pValues() { return statistics.pValues(); }
    public double[] pValues(PValueScale scale) {
        return statistics.pValues(scale);
    }
    public double[] log10PValues() { return statistics.log10PValues(); }
    public double[] negativeLog10PValues() {
        return statistics.negativeLog10PValues();
    }
    public RemlResult nullModel() { return nullModel; }
}
