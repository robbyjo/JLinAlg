/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.reml.RemlResult;

/** Gaussian GAMM result with smooth and covariance-component contributions. */
public final class GammResult {
    private final RemlResult reml;
    private final int parametricColumns;
    private final List<SmoothTermEstimate> smoothTerms;
    private final Map<String, double[]> randomContributions;
    private final double[] fittedValues;
    private final double[] residuals;
    private final double totalEffectiveDegreesOfFreedom;

    GammResult(
            RemlResult reml,
            int parametricColumns,
            List<SmoothTermEstimate> smoothTerms,
            Map<String, double[]> randomContributions,
            double[] fittedValues,
            double[] residuals,
            double totalEffectiveDegreesOfFreedom) {
        this.reml = reml;
        this.parametricColumns = parametricColumns;
        this.smoothTerms = List.copyOf(smoothTerms);
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomContributions.forEach(
            (name, values) -> copied.put(name, values.clone()));
        this.randomContributions = java.util.Collections.unmodifiableMap(copied);
        this.fittedValues = fittedValues.clone();
        this.residuals = residuals.clone();
        this.totalEffectiveDegreesOfFreedom = totalEffectiveDegreesOfFreedom;
    }

    public RemlResult reml() { return reml; }
    public double[] parametricCoefficients() {
        return prefix(reml.beta(), parametricColumns);
    }
    public double[] parametricStandardErrors() {
        return prefix(reml.standardErrors(), parametricColumns);
    }
    public double[] parametricTStatistics() {
        return prefix(reml.tStatistics(), parametricColumns);
    }
    public double[] parametricPValues() {
        return prefix(reml.pValues(), parametricColumns);
    }
    public List<SmoothTermEstimate> smoothTerms() { return smoothTerms; }
    public Map<String, double[]> randomContributions() {
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomContributions.forEach(
            (name, values) -> copied.put(name, values.clone()));
        return java.util.Collections.unmodifiableMap(copied);
    }
    public double[] randomContribution(String name) {
        double[] result = randomContributions.get(name);
        if (result == null) {
            throw new IllegalArgumentException(
                "unknown covariance component: " + name);
        }
        return result.clone();
    }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
    public double totalEffectiveDegreesOfFreedom() {
        return totalEffectiveDegreesOfFreedom;
    }
    public double residualDegreesOfFreedom() {
        return reml.observations() - totalEffectiveDegreesOfFreedom;
    }

    private static double[] prefix(double[] values, int length) {
        return java.util.Arrays.copyOf(values, length);
    }
}
