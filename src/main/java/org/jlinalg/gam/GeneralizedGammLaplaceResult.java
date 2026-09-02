/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.LinkedHashMap;
import java.util.Map;
import org.jlinalg.glmm.GlmmLaplaceResult;

/** Laplace GAMM result separating smooth and non-smooth covariance contributions. */
public final class GeneralizedGammLaplaceResult {
    private final GlmmLaplaceResult fit;
    private final int parametricColumns;
    private final Map<String, double[]> smoothContributions;
    private final Map<String, Double> smoothingParameters;
    private final Map<String, double[]> randomContributions;

    GeneralizedGammLaplaceResult(
            GlmmLaplaceResult fit,
            int parametricColumns,
            Map<String, double[]> smoothContributions,
            Map<String, Double> smoothingParameters,
            Map<String, double[]> randomContributions) {
        this.fit = fit;
        this.parametricColumns = parametricColumns;
        this.smoothContributions = copy(smoothContributions);
        this.smoothingParameters = Map.copyOf(smoothingParameters);
        this.randomContributions = copy(randomContributions);
    }
    public GlmmLaplaceResult mixedModel() { return fit; }
    public double[] parametricCoefficients() {
        return java.util.Arrays.copyOf(fit.beta(), parametricColumns);
    }
    public double[] parametricStandardErrors() {
        return java.util.Arrays.copyOf(fit.standardErrors(), parametricColumns);
    }
    public double[] parametricStatistics() {
        return java.util.Arrays.copyOf(fit.statistics(), parametricColumns);
    }
    public double[] parametricPValues() {
        return java.util.Arrays.copyOf(fit.pValues(), parametricColumns);
    }
    public Map<String, double[]> smoothContributions() {
        return copy(smoothContributions);
    }
    public Map<String, Double> smoothingParameters() {
        return smoothingParameters;
    }
    public Map<String, double[]> randomContributions() {
        return copy(randomContributions);
    }
    public double[] fittedMeans() { return fit.fittedMeans(); }
    public double marginalLogLikelihood() { return fit.marginalLogLikelihood(); }
    public boolean converged() { return fit.converged(); }
    public String convergenceMessage() { return fit.convergenceMessage(); }

    private static Map<String, double[]> copy(Map<String, double[]> source) {
        Map<String, double[]> result = new LinkedHashMap<>();
        source.forEach((name, values) -> result.put(name, values.clone()));
        return java.util.Collections.unmodifiableMap(result);
    }
}
