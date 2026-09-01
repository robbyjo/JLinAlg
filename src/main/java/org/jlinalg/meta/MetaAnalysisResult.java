/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Pooled effect and heterogeneity diagnostics from a univariate meta-analysis. */
public record MetaAnalysisResult(
        List<String> studyNames,
        MetaAnalysisMethod method,
        TauSquaredEstimator tauSquaredEstimator,
        MetaInferenceMethod inferenceMethod,
        AssociationStatistics associationStatistics,
        double confidenceLower,
        double confidenceUpper,
        double predictionLower,
        double predictionUpper,
        double cochranQ,
        double cochranQDegreesOfFreedom,
        double cochranQPValue,
        double tauSquared,
        double iSquared,
        double hSquared,
        double[] normalizedWeights,
        BackendProvenance backend) {

    public MetaAnalysisResult {
        studyNames = List.copyOf(studyNames);
        normalizedWeights = normalizedWeights.clone();
    }

    public double pooledEffectSize() { return associationStatistics.beta()[0]; }
    public double standardError() { return associationStatistics.standardErrors()[0]; }
    public double statistic() { return associationStatistics.statistics()[0]; }
    public double pValue() { return associationStatistics.pValues()[0]; }
    public double log10PValue() { return associationStatistics.log10PValues()[0]; }
    public double negativeLog10PValue() {
        return associationStatistics.negativeLog10PValues()[0];
    }
    public double tau() { return Math.sqrt(tauSquared); }

    @Override public double[] normalizedWeights() { return normalizedWeights.clone(); }
}
