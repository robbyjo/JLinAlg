/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.glmm.GlmmPqlResult;

/** Generalized GAMM result with smooth and covariance-component modes. */
public final class GeneralizedGammResult {
    private final GlmmPqlResult workingModel;
    private final int parametricColumns;
    private final List<SmoothTermEstimate> smoothTerms;
    private final Map<String, double[]> randomLinearPredictors;
    private final double totalEffectiveDegreesOfFreedom;

    GeneralizedGammResult(
            GlmmPqlResult workingModel,
            int parametricColumns,
            List<SmoothTermEstimate> smoothTerms,
            Map<String, double[]> randomLinearPredictors,
            double totalEffectiveDegreesOfFreedom) {
        this.workingModel = workingModel;
        this.parametricColumns = parametricColumns;
        this.smoothTerms = List.copyOf(smoothTerms);
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomLinearPredictors.forEach(
            (name, values) -> copied.put(name, values.clone()));
        this.randomLinearPredictors =
            java.util.Collections.unmodifiableMap(copied);
        this.totalEffectiveDegreesOfFreedom = totalEffectiveDegreesOfFreedom;
    }

    public GlmmPqlResult workingModel() { return workingModel; }
    public boolean converged() { return workingModel.converged(); }
    public String convergenceMessage() {
        return workingModel.convergenceMessage();
    }
    public double[] parametricCoefficients() {
        return prefix(workingModel.beta(), parametricColumns);
    }
    public double[] parametricStandardErrors() {
        return prefix(workingModel.standardErrors(), parametricColumns);
    }
    public double[] parametricTStatistics() {
        return prefix(workingModel.tStatistics(), parametricColumns);
    }
    public double[] parametricPValues() {
        return prefix(workingModel.pValues(), parametricColumns);
    }
    public double[] fittedMeans() { return workingModel.fittedMeans(); }
    public double[] linearPredictor() {
        return workingModel.linearPredictor();
    }
    public List<SmoothTermEstimate> smoothTerms() { return smoothTerms; }
    public double[] randomLinearPredictor(String name) {
        double[] result = randomLinearPredictors.get(name);
        if (result == null) {
            throw new IllegalArgumentException(
                "unknown covariance component: " + name);
        }
        return result.clone();
    }
    public Map<String, double[]> randomLinearPredictors() {
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomLinearPredictors.forEach(
            (name, values) -> copied.put(name, values.clone()));
        return java.util.Collections.unmodifiableMap(copied);
    }
    public double totalEffectiveDegreesOfFreedom() {
        return totalEffectiveDegreesOfFreedom;
    }

    private static double[] prefix(double[] values, int length) {
        return java.util.Arrays.copyOf(values, length);
    }
}
