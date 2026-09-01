/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.reml.VarianceEstimation;

/** Sparse-equation Gaussian mixed-model estimates and diagnostics. */
public final class SparseLinearMixedModelResult {
    private final List<String> componentNames;
    private final double[] varianceComponents;
    private final AssociationStatistics associationStatistics;
    private final double[] fixedEffectCovariance;
    private final List<RandomEffectEstimates> randomEffects;
    private final Map<String, RandomEffectEstimates> randomByName;
    private final double[] conditionalFittedValues;
    private final double[] conditionalResiduals;
    private final double logLikelihood;
    private final VarianceEstimation varianceEstimation;
    private final int functionEvaluations;
    private final boolean converged;
    private final int randomCoefficientCount;
    private final int equationNonzeroCount;
    private final int factorNonzeroCount;
    private final BackendProvenance backend;

    SparseLinearMixedModelResult(
            List<String> componentNames,
            double[] varianceComponents,
            AssociationStatistics associationStatistics,
            double[] fixedEffectCovariance,
            List<RandomEffectEstimates> randomEffects,
            double[] conditionalFittedValues,
            double[] conditionalResiduals,
            double logLikelihood,
            VarianceEstimation varianceEstimation,
            int functionEvaluations,
            boolean converged,
            int randomCoefficientCount,
            int equationNonzeroCount,
            int factorNonzeroCount,
            BackendProvenance backend) {
        this.componentNames = List.copyOf(componentNames);
        this.varianceComponents = varianceComponents.clone();
        this.associationStatistics = associationStatistics;
        this.fixedEffectCovariance = fixedEffectCovariance.clone();
        this.randomEffects = List.copyOf(randomEffects);
        Map<String, RandomEffectEstimates> indexed = new LinkedHashMap<>();
        for (RandomEffectEstimates value : randomEffects)
            indexed.put(value.termName(), value);
        this.randomByName = Map.copyOf(indexed);
        this.conditionalFittedValues = conditionalFittedValues.clone();
        this.conditionalResiduals = conditionalResiduals.clone();
        this.logLikelihood = logLikelihood;
        this.varianceEstimation = varianceEstimation;
        this.functionEvaluations = functionEvaluations;
        this.converged = converged;
        this.randomCoefficientCount = randomCoefficientCount;
        this.equationNonzeroCount = equationNonzeroCount;
        this.factorNonzeroCount = factorNonzeroCount;
        this.backend = backend;
    }

    public List<String> componentNames() { return componentNames; }
    public double[] varianceComponents() { return varianceComponents.clone(); }
    public AssociationStatistics associationStatistics() {
        return associationStatistics;
    }
    public double[] beta() { return associationStatistics.beta(); }
    public double[] fixef() { return associationStatistics.beta(); }
    public double[] effectSizes() { return associationStatistics.effectSizes(); }
    public double[] standardErrors() {
        return associationStatistics.standardErrors();
    }
    public double[] tStatistics() {
        return associationStatistics.statistics();
    }
    public double[] pValues() { return associationStatistics.pValues(); }
    public double[] fixedEffectCovariance() {
        return fixedEffectCovariance.clone();
    }
    public List<RandomEffectEstimates> randomEffects() { return randomEffects; }
    public Map<String, RandomEffectEstimates> ranef() { return randomByName; }
    public List<VarianceComponentSummary> varCorr() {
        java.util.ArrayList<VarianceComponentSummary> result =
            new java.util.ArrayList<>(varianceComponents.length);
        for (int index = 0; index < varianceComponents.length; index++)
            result.add(VarianceComponentSummary.of(
                componentNames.get(index), varianceComponents[index]));
        return List.copyOf(result);
    }
    public RandomEffectEstimates randomEffects(String termName) {
        RandomEffectEstimates result = randomByName.get(termName);
        if (result == null)
            throw new IllegalArgumentException(
                "unknown random-effect term: " + termName);
        return result;
    }
    public double[] conditionalFittedValues() {
        return conditionalFittedValues.clone();
    }
    public double[] fittedValues() { return conditionalFittedValues(); }
    public double[] conditionalResiduals() {
        return conditionalResiduals.clone();
    }
    public double[] residuals() { return conditionalResiduals(); }
    public double logLikelihood() { return logLikelihood; }
    public VarianceEstimation varianceEstimation() { return varianceEstimation; }
    public int functionEvaluations() { return functionEvaluations; }
    public boolean converged() { return converged; }
    public int randomCoefficientCount() { return randomCoefficientCount; }
    public int equationNonzeroCount() { return equationNonzeroCount; }
    public int factorNonzeroCount() { return factorNonzeroCount; }
    public BackendProvenance backend() { return backend; }
}
