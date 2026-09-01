/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.Normal;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Laplace-approximated Gaussian frailty Cox model result. */
public final class CoxMixedResult {
    private final double[] beta;
    private final double[] covariance;
    private final AssociationStatistics statistics;
    private final double[] hazardRatioLower;
    private final double[] hazardRatioUpper;
    private final List<CoxRandomEffectEstimates> randomEffects;
    private final Map<String, CoxRandomEffectEstimates> randomEffectsByName;
    private final List<BaselineHazardPoint> baselineHazard;
    private final double partialLogLikelihood;
    private final double penalizedLogLikelihood;
    private final double laplaceLogLikelihood;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final BackendProvenance backend;

    CoxMixedResult(
            double[] beta, double[] covariance,
            List<CoxRandomEffectEstimates> randomEffects,
            List<BaselineHazardPoint> baselineHazard,
            double partialLogLikelihood, double penalizedLogLikelihood,
            double laplaceLogLikelihood, CoxMixedOptions options,
            int iterations, boolean converged, String convergenceMessage,
            BackendProvenance backend) {
        this.beta = beta.clone();
        this.covariance = covariance.clone();
        double[] standardErrors = new double[beta.length];
        for (int column = 0; column < beta.length; column++)
            standardErrors[column] = Math.sqrt(Math.max(0,
                covariance[column * beta.length + column]));
        statistics = AssociationStatistics.normal(beta, standardErrors);
        double critical = Normal.quantile(
            0.5 + options.coxOptions().confidenceLevel() / 2,
            0, 1, true, false);
        hazardRatioLower = new double[beta.length];
        hazardRatioUpper = new double[beta.length];
        for (int column = 0; column < beta.length; column++) {
            hazardRatioLower[column] = Math.exp(beta[column]
                - critical * standardErrors[column]);
            hazardRatioUpper[column] = Math.exp(beta[column]
                + critical * standardErrors[column]);
        }
        this.randomEffects = List.copyOf(randomEffects);
        Map<String, CoxRandomEffectEstimates> byName = new LinkedHashMap<>();
        for (CoxRandomEffectEstimates value : randomEffects)
            if (byName.put(value.termName(), value) != null)
                throw new IllegalArgumentException(
                    "duplicate Cox frailty term: " + value.termName());
        randomEffectsByName = Map.copyOf(byName);
        this.baselineHazard = List.copyOf(baselineHazard);
        this.partialLogLikelihood = partialLogLikelihood;
        this.penalizedLogLikelihood = penalizedLogLikelihood;
        this.laplaceLogLikelihood = laplaceLogLikelihood;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
        this.backend = backend;
    }

    public double[] beta() { return beta.clone(); }
    public double[] fixef() { return beta(); }
    public double[] hazardRatios() {
        return java.util.Arrays.stream(beta).map(Math::exp).toArray();
    }
    public double[] standardErrors() { return statistics.standardErrors(); }
    public double[] zStatistics() { return statistics.statistics(); }
    public double[] pValues() { return statistics.pValues(); }
    public double[] negativeLog10PValues() {
        return statistics.negativeLog10PValues();
    }
    public double[] hazardRatioConfidenceLower() {
        return hazardRatioLower.clone();
    }
    public double[] hazardRatioConfidenceUpper() {
        return hazardRatioUpper.clone();
    }
    public double[] covariance() { return covariance.clone(); }
    public AssociationStatistics associationStatistics() { return statistics; }
    public List<CoxRandomEffectEstimates> randomEffects() {
        return randomEffects;
    }
    public Map<String, CoxRandomEffectEstimates> ranef() {
        return randomEffectsByName;
    }
    public CoxRandomEffectEstimates randomEffects(String name) {
        CoxRandomEffectEstimates value = randomEffectsByName.get(name);
        if (value == null)
            throw new IllegalArgumentException("unknown Cox frailty term: " + name);
        return value;
    }
    public List<BaselineHazardPoint> baselineHazard() { return baselineHazard; }
    public double partialLogLikelihood() { return partialLogLikelihood; }
    public double penalizedLogLikelihood() { return penalizedLogLikelihood; }
    public double laplaceLogLikelihood() { return laplaceLogLikelihood; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public BackendProvenance backend() { return backend; }
}
