/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/**
 * Columnar results from a prepared batch of independent meta-analyses.
 * Array positions follow the analysis-row order supplied at preparation time.
 */
public final class MetaAnalysisBatchResult {
    private final MetaAnalysisMethod method;
    private final TauSquaredEstimator tauSquaredEstimator;
    private final MetaInferenceMethod inferenceMethod;
    private final double[] pooledEffectSizes;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] negativeLog10PValues;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final double[] predictionLower;
    private final double[] predictionUpper;
    private final double[] cochranQ;
    private final double[] cochranQPValues;
    private final double[] tauSquared;
    private final double[] iSquared;
    private final double[] hSquared;

    MetaAnalysisBatchResult(
            MetaAnalysisMethod method,
            TauSquaredEstimator tauSquaredEstimator,
            MetaInferenceMethod inferenceMethod,
            double[] pooledEffectSizes,
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] negativeLog10PValues,
            double[] confidenceLower,
            double[] confidenceUpper,
            double[] predictionLower,
            double[] predictionUpper,
            double[] cochranQ,
            double[] cochranQPValues,
            double[] tauSquared,
            double[] iSquared,
            double[] hSquared) {
        this.method = method;
        this.tauSquaredEstimator = tauSquaredEstimator;
        this.inferenceMethod = inferenceMethod;
        this.pooledEffectSizes = pooledEffectSizes;
        this.standardErrors = standardErrors;
        this.statistics = statistics;
        this.pValues = pValues;
        this.negativeLog10PValues = negativeLog10PValues;
        this.confidenceLower = confidenceLower;
        this.confidenceUpper = confidenceUpper;
        this.predictionLower = predictionLower;
        this.predictionUpper = predictionUpper;
        this.cochranQ = cochranQ;
        this.cochranQPValues = cochranQPValues;
        this.tauSquared = tauSquared;
        this.iSquared = iSquared;
        this.hSquared = hSquared;
    }

    public MetaAnalysisMethod method() { return method; }
    public TauSquaredEstimator tauSquaredEstimator() { return tauSquaredEstimator; }
    public MetaInferenceMethod inferenceMethod() { return inferenceMethod; }
    public int analyses() { return pooledEffectSizes.length; }
    public double[] pooledEffectSizes() { return pooledEffectSizes.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] negativeLog10PValues() {
        return negativeLog10PValues.clone();
    }
    public double[] log10PValues() {
        double[] result = negativeLog10PValues.clone();
        for (int index = 0; index < result.length; index++)
            result[index] = -result[index];
        return result;
    }
    public double[] confidenceLower() { return confidenceLower.clone(); }
    public double[] confidenceUpper() { return confidenceUpper.clone(); }
    public double[] predictionLower() { return predictionLower.clone(); }
    public double[] predictionUpper() { return predictionUpper.clone(); }
    public double[] cochranQ() { return cochranQ.clone(); }
    public double[] cochranQPValues() { return cochranQPValues.clone(); }
    public double[] tauSquared() { return tauSquared.clone(); }
    public double[] iSquared() { return iSquared.clone(); }
    public double[] hSquared() { return hSquared.clone(); }
}
