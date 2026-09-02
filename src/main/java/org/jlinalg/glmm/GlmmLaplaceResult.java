/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.inference.AssociationStatistics;

/** Conditional modes and marginal likelihood from a Laplace GLMM fit. */
public final class GlmmLaplaceResult {
    private final String family;
    private final List<String> componentNames;
    private final double[] variances;
    private final AssociationStatistics association;
    private final double[] fixedCovariance;
    private final Map<String, double[]> componentPredictors;
    private final double[] linearPredictor;
    private final double[] fittedMeans;
    private final double logLikelihood;
    private final int outerIterations;
    private final int modeIterations;
    private final boolean converged;

    GlmmLaplaceResult(
            String family,
            List<String> componentNames,
            double[] variances,
            AssociationStatistics association,
            double[] fixedCovariance,
            Map<String, double[]> componentPredictors,
            double[] linearPredictor,
            double[] fittedMeans,
            double logLikelihood,
            int outerIterations,
            int modeIterations,
            boolean converged) {
        this.family = family;
        this.componentNames = List.copyOf(componentNames);
        this.variances = variances.clone();
        this.association = association;
        this.fixedCovariance = fixedCovariance.clone();
        Map<String, double[]> copied = new LinkedHashMap<>();
        componentPredictors.forEach((name, values) ->
            copied.put(name, values.clone()));
        this.componentPredictors = java.util.Collections.unmodifiableMap(copied);
        this.linearPredictor = linearPredictor.clone();
        this.fittedMeans = fittedMeans.clone();
        this.logLikelihood = logLikelihood;
        this.outerIterations = outerIterations;
        this.modeIterations = modeIterations;
        this.converged = converged;
    }
    public String family() { return family; }
    public List<String> componentNames() { return componentNames; }
    public double[] varianceComponents() { return variances.clone(); }
    public AssociationStatistics associationStatistics() { return association; }
    public double[] beta() { return association.beta(); }
    public double[] standardErrors() { return association.standardErrors(); }
    public double[] statistics() { return association.statistics(); }
    public double[] pValues() { return association.pValues(); }
    public double[] fixedEffectCovariance() { return fixedCovariance.clone(); }
    public double[] componentPredictor(String name) {
        double[] values = componentPredictors.get(name);
        if (values == null) throw new IllegalArgumentException("unknown component: " + name);
        return values.clone();
    }
    public Map<String, double[]> componentPredictors() {
        Map<String, double[]> result = new LinkedHashMap<>();
        componentPredictors.forEach((name, values) ->
            result.put(name, values.clone()));
        return java.util.Collections.unmodifiableMap(result);
    }
    public double[] linearPredictor() { return linearPredictor.clone(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double marginalLogLikelihood() { return logLikelihood; }
    public int outerIterations() { return outerIterations; }
    public int modeIterations() { return modeIterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() {
        return converged ? "Laplace variance and conditional-mode tolerances reached"
            : "Laplace optimizer stopped before all tolerances were reached";
    }
}
