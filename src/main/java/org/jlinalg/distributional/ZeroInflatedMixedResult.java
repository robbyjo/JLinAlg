/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Maximum-Laplace estimates for a zero-inflated count mixed model. */
public final class ZeroInflatedMixedResult {
    private final String family;
    private final double[] countCoefficients;
    private final double[] zeroCoefficients;
    private final double[] dispersionCoefficients;
    private final List<String> componentNames;
    private final double[] variances;
    private final Map<String, Double> correlations;
    private final Map<String, double[]> randomEffects;
    private final double[] conditionalCountMeans;
    private final double[] structuralZeroProbabilities;
    private final double[] sizes;
    private final double[] fittedMeans;
    private final double[] fittedZeroProbabilities;
    private final double marginalLogLikelihood;
    private final int objectiveEvaluations;
    private final int modeIterations;
    private final boolean converged;
    private final int randomCoefficientCount;
    private final int sparseEquationNonzeroCount;
    private final int factorNonzeroCount;
    private final List<String> marginalParameterNames;
    private final double[] marginalCovariance;
    private final List<String> warnings;
    private final List<String> outerParameterNames;
    private final double[] outerParameterEstimates;

    ZeroInflatedMixedResult(
            String family,
            double[] countCoefficients,
            double[] zeroCoefficients,
            double[] dispersionCoefficients,
            List<String> componentNames,
            double[] variances,
            Map<String, Double> correlations,
            Map<String, double[]> randomEffects,
            double[] conditionalCountMeans,
            double[] structuralZeroProbabilities,
            double[] sizes,
            double[] fittedMeans,
            double[] fittedZeroProbabilities,
            double marginalLogLikelihood,
            int objectiveEvaluations,
            int modeIterations,
            boolean converged,
            int randomCoefficientCount,
            int sparseEquationNonzeroCount,
            int factorNonzeroCount,
            List<String> marginalParameterNames,
            double[] marginalCovariance,
            List<String> warnings,
            List<String> outerParameterNames,
            double[] outerParameterEstimates) {
        this.family = family;
        this.countCoefficients = countCoefficients.clone();
        this.zeroCoefficients = zeroCoefficients.clone();
        this.dispersionCoefficients = dispersionCoefficients.clone();
        this.componentNames = List.copyOf(componentNames);
        this.variances = variances.clone();
        this.correlations = Map.copyOf(correlations);
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomEffects.forEach((name, values) ->
            copied.put(name, values.clone()));
        this.randomEffects = java.util.Collections.unmodifiableMap(copied);
        this.conditionalCountMeans = conditionalCountMeans.clone();
        this.structuralZeroProbabilities =
            structuralZeroProbabilities.clone();
        this.sizes = sizes.clone();
        this.fittedMeans = fittedMeans.clone();
        this.fittedZeroProbabilities = fittedZeroProbabilities.clone();
        this.marginalLogLikelihood = marginalLogLikelihood;
        this.objectiveEvaluations = objectiveEvaluations;
        this.modeIterations = modeIterations;
        this.converged = converged;
        this.randomCoefficientCount = randomCoefficientCount;
        this.sparseEquationNonzeroCount = sparseEquationNonzeroCount;
        this.factorNonzeroCount = factorNonzeroCount;
        this.marginalParameterNames = List.copyOf(marginalParameterNames);
        this.marginalCovariance = marginalCovariance.clone();
        this.warnings = List.copyOf(warnings);
        this.outerParameterNames = List.copyOf(outerParameterNames);
        this.outerParameterEstimates = outerParameterEstimates.clone();
    }

    public String family() { return family; }
    public double[] countCoefficients() { return countCoefficients.clone(); }
    public double[] zeroCoefficients() { return zeroCoefficients.clone(); }
    public double[] dispersionCoefficients() {
        return dispersionCoefficients.clone();
    }
    public List<String> componentNames() { return componentNames; }
    public double[] varianceComponents() { return variances.clone(); }
    /** Estimated count/zero correlations for named two-process effects. */
    public Map<String, Double> correlations() { return correlations; }
    public Map<String, double[]> randomEffects() {
        Map<String, double[]> result = new LinkedHashMap<>();
        randomEffects.forEach((name, values) ->
            result.put(name, values.clone()));
        return java.util.Collections.unmodifiableMap(result);
    }
    public double[] randomEffects(String component) {
        double[] result = randomEffects.get(component);
        if (result == null) {
            throw new IllegalArgumentException(
                "unknown random-effect component: " + component);
        }
        return result.clone();
    }
    public double[] conditionalCountMeans() {
        return conditionalCountMeans.clone();
    }
    public double[] structuralZeroProbabilities() {
        return structuralZeroProbabilities.clone();
    }
    /** NB2 size per observation; empty for a Poisson model. */
    public double[] sizes() { return sizes.clone(); }
    /** Marginal response means conditional on fitted random-effect modes. */
    public double[] fittedMeans() { return fittedMeans.clone(); }
    /** Total fitted zero mass, including sampling and structural zeros. */
    public double[] fittedZeroProbabilities() {
        return fittedZeroProbabilities.clone();
    }
    public double marginalLogLikelihood() { return marginalLogLikelihood; }
    public int objectiveEvaluations() { return objectiveEvaluations; }
    public int modeIterations() { return modeIterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() {
        return converged
            ? "outer optimizer and random-effect mode tolerances reached"
            : "outer optimizer or random-effect mode stopped before convergence";
    }
    public int randomCoefficientCount() { return randomCoefficientCount; }
    public int sparseEquationNonzeroCount() {
        return sparseEquationNonzeroCount;
    }
    public int factorNonzeroCount() { return factorNonzeroCount; }
    /** Names for the count, zero, and log-size marginal covariance. */
    public List<String> marginalParameterNames() {
        return marginalParameterNames;
    }
    /** Row-major observed-Hessian covariance; empty unless requested. */
    public double[] marginalCovariance() { return marginalCovariance.clone(); }
    public boolean inferenceAvailable() { return marginalCovariance.length > 0; }
    public double[] standardErrors() {
        int dimension = marginalParameterNames.size();
        double[] result = new double[dimension];
        for (int index = 0; index < dimension; index++) {
            double variance = marginalCovariance[index * dimension + index];
            result[index] = variance >= 0.0 ? Math.sqrt(variance) : Double.NaN;
        }
        return result;
    }
    /** Convergence, singularity, and boundary diagnostics. */
    public List<String> warnings() { return warnings; }
    public boolean singular() {
        return warnings.stream().anyMatch(value ->
            value.contains("boundary") || value.contains("singular"));
    }
    /** Optimizer-scale names, including log variances and Fisher-z terms. */
    public List<String> outerParameterNames() { return outerParameterNames; }
    /** Optimizer-scale estimates in {@link #outerParameterNames()} order. */
    public double[] outerParameterEstimates() {
        return outerParameterEstimates.clone();
    }
}
