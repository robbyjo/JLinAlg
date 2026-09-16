/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;
import org.jlinalg.inference.AssociationStatistics;

/** Covariance-aware direct effects for multiple exposures and correlated outcomes. */
public final class MultivariateMrResult {
    private final List<String> exposureNames;
    private final List<String> outcomeNames;
    private final AssociationStatistics statistics;
    private final double[] covariance;
    private final double[] marginalFStatistics;
    private final MultivariateMrJointTest overallTest;
    private final List<MultivariateMrJointTest> exposureTests;
    private final List<MultivariateMrJointTest> outcomeTests;
    private final double heterogeneity;
    private final int heterogeneityDegreesOfFreedom;
    private final double heterogeneityPValue;
    private final double dispersion;
    private final boolean multiplicativeRandomEffects;

    MultivariateMrResult(List<String> exposureNames, List<String> outcomeNames,
            AssociationStatistics statistics, double[] covariance,
            double[] marginalFStatistics, MultivariateMrJointTest overallTest,
            List<MultivariateMrJointTest> exposureTests,
            List<MultivariateMrJointTest> outcomeTests,
            double heterogeneity, int heterogeneityDegreesOfFreedom,
            double heterogeneityPValue, double dispersion,
            boolean multiplicativeRandomEffects) {
        this.exposureNames = List.copyOf(exposureNames);
        this.outcomeNames = List.copyOf(outcomeNames);
        this.statistics = statistics;
        this.covariance = covariance.clone();
        this.marginalFStatistics = marginalFStatistics.clone();
        this.overallTest = overallTest;
        this.exposureTests = List.copyOf(exposureTests);
        this.outcomeTests = List.copyOf(outcomeTests);
        this.heterogeneity = heterogeneity;
        this.heterogeneityDegreesOfFreedom = heterogeneityDegreesOfFreedom;
        this.heterogeneityPValue = heterogeneityPValue;
        this.dispersion = dispersion;
        this.multiplicativeRandomEffects = multiplicativeRandomEffects;
    }

    public List<String> exposureNames() { return exposureNames; }
    public List<String> outcomeNames() { return outcomeNames; }
    public AssociationStatistics statistics() { return statistics; }
    /** Coefficients in outcome-major order, with exposures varying fastest. */
    public double[] beta() { return statistics.beta(); }
    public double[] standardErrors() { return statistics.standardErrors(); }
    public double[] pValues() { return statistics.pValues(); }
    public double[] covariance() { return covariance.clone(); }
    public double[] marginalFStatistics() { return marginalFStatistics.clone(); }
    public MultivariateMrJointTest overallTest() { return overallTest; }
    public List<MultivariateMrJointTest> exposureJointTests() { return exposureTests; }
    public List<MultivariateMrJointTest> outcomeJointTests() { return outcomeTests; }
    public double cochranQ() { return heterogeneity; }
    public int heterogeneityDegreesOfFreedom() { return heterogeneityDegreesOfFreedom; }
    public double heterogeneityPValue() { return heterogeneityPValue; }
    public double dispersion() { return dispersion; }
    public boolean multiplicativeRandomEffects() { return multiplicativeRandomEffects; }

    public int coefficientIndex(int exposure, int outcome) {
        if (exposure < 0 || exposure >= exposureNames.size()
                || outcome < 0 || outcome >= outcomeNames.size())
            throw new IndexOutOfBoundsException("exposure or outcome index is invalid");
        return outcome * exposureNames.size() + exposure;
    }

    public double effect(int exposure, int outcome) {
        return beta()[coefficientIndex(exposure, outcome)];
    }
}
