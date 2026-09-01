/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Coefficients, moderator test, and residual heterogeneity from meta-regression. */
public record MetaRegressionResult(
        List<String> coefficientNames,
        MetaAnalysisMethod method,
        TauSquaredEstimator tauSquaredEstimator,
        MetaInferenceMethod inferenceMethod,
        AssociationStatistics associationStatistics,
        double[] coefficientCovariance,
        double residualQ,
        double residualQDegreesOfFreedom,
        double residualQPValue,
        double moderatorQ,
        double moderatorQDegreesOfFreedom,
        double moderatorQPValue,
        double tauSquared,
        double residualISquared,
        double residualHSquared,
        double heterogeneityRSquared,
        double[] normalizedWeights,
        BackendProvenance backend) {

    public MetaRegressionResult {
        coefficientNames = List.copyOf(coefficientNames);
        coefficientCovariance = coefficientCovariance.clone();
        normalizedWeights = normalizedWeights.clone();
    }

    public double[] beta() { return associationStatistics.beta(); }
    public double[] effectSizes() { return associationStatistics.effectSizes(); }
    public double[] standardErrors() { return associationStatistics.standardErrors(); }
    public double[] statistics() { return associationStatistics.statistics(); }
    public double[] pValues() { return associationStatistics.pValues(); }
    public double[] log10PValues() { return associationStatistics.log10PValues(); }
    public double[] negativeLog10PValues() {
        return associationStatistics.negativeLog10PValues();
    }
    public double tau() { return Math.sqrt(tauSquared); }

    @Override public double[] coefficientCovariance() {
        return coefficientCovariance.clone();
    }
    @Override public double[] normalizedWeights() { return normalizedWeights.clone(); }
}
