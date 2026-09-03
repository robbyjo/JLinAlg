/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.ContrastTestResult;
import org.jlinalg.inference.LinearHypothesis;

/** Immutable estimates and convergence metadata from a Gaussian REML fit. */
public final class RemlResult {
    private final List<String> componentNames;
    private final double[] varianceComponents;
    private final double[] fixedEffects;
    private final double[] fixedEffectCovariance;
    private final double[] fixedEffectInferenceCovariance;
    private final double[] fixedEffectStandardErrors;
    private final double[] fittedValues;
    private final double[] residuals;
    private final double[] logVarianceScore;
    private final AssociationStatistics associationStatistics;
    private final double restrictedLogLikelihood;
    private final VarianceEstimation varianceEstimation;
    private final int observations;
    private final int fixedEffectCount;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final BackendProvenance backend;

    RemlResult(
            List<String> componentNames,
            double[] varianceComponents,
            double[] fixedEffects,
            double[] fixedEffectCovariance,
            double[] fixedEffectInferenceCovariance,
            double[] fixedEffectStandardErrors,
            double[] fittedValues,
            double[] residuals,
            double[] logVarianceScore,
            AssociationStatistics associationStatistics,
            double restrictedLogLikelihood,
            VarianceEstimation varianceEstimation,
            int observations,
            int fixedEffectCount,
            int iterations,
            boolean converged,
            String convergenceMessage,
            BackendProvenance backend) {
        this.componentNames = List.copyOf(componentNames);
        this.varianceComponents = varianceComponents.clone();
        this.fixedEffects = fixedEffects.clone();
        this.fixedEffectCovariance = fixedEffectCovariance.clone();
        this.fixedEffectInferenceCovariance =
            fixedEffectInferenceCovariance.clone();
        this.fixedEffectStandardErrors = fixedEffectStandardErrors.clone();
        this.fittedValues = fittedValues.clone();
        this.residuals = residuals.clone();
        this.logVarianceScore = logVarianceScore.clone();
        this.associationStatistics = Objects.requireNonNull(
            associationStatistics, "associationStatistics");
        this.restrictedLogLikelihood = restrictedLogLikelihood;
        this.varianceEstimation = Objects.requireNonNull(
            varianceEstimation, "varianceEstimation");
        this.observations = observations;
        this.fixedEffectCount = fixedEffectCount;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = Objects.requireNonNull(
            convergenceMessage, "convergenceMessage");
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    /**
     * Adapts an exact coefficient-space Gaussian fit to the common REML
     * result contract. Score values are unavailable from derivative-free
     * coefficient-space optimization and are reported as zero.
     */
    public static RemlResult fromCoefficientSpace(
            List<String> componentNames,
            double[] varianceComponents,
            AssociationStatistics associationStatistics,
            double[] fixedEffectCovariance,
            double[] fittedValues,
            double[] residuals,
            double logLikelihood,
            VarianceEstimation varianceEstimation,
            int observations,
            int fixedEffectCount,
            int functionEvaluations,
            boolean converged,
            BackendProvenance backend) {
        Objects.requireNonNull(
            associationStatistics, "associationStatistics");
        return new RemlResult(componentNames, varianceComponents,
            associationStatistics.beta(), fixedEffectCovariance,
            fixedEffectCovariance, associationStatistics.standardErrors(),
            fittedValues, residuals, new double[varianceComponents.length],
            associationStatistics, logLikelihood, varianceEstimation,
            observations, fixedEffectCount, functionEvaluations, converged,
            converged
                ? "coefficient-space optimizer converged"
                : "coefficient-space optimizer reached its evaluation limit",
            backend);
    }

    public List<String> componentNames() { return componentNames; }
    public double[] varianceComponents() { return varianceComponents.clone(); }
    public double[] fixedEffects() { return fixedEffects.clone(); }
    /** Model-based covariance treating fitted variance components as fixed. */
    public double[] fixedEffectCovariance() { return fixedEffectCovariance.clone(); }
    /** Covariance used for association inference; KR adjusts this matrix. */
    public double[] fixedEffectInferenceCovariance() {
        return fixedEffectInferenceCovariance.clone();
    }
    public double[] fixedEffectStandardErrors() { return fixedEffectStandardErrors.clone(); }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
    public double[] logVarianceScore() { return logVarianceScore.clone(); }
    /** Returns beta, SE, t statistics, denominator DF, and two-sided p-values. */
    public AssociationStatistics associationStatistics() {
        return associationStatistics;
    }
    /** Alias for fixed effects, convenient for association scans. */
    public double[] beta() { return associationStatistics.beta(); }
    public double[] standardErrors() {
        return associationStatistics.standardErrors();
    }
    public double[] statistics() { return associationStatistics.statistics(); }
    public double[] tStatistics() { return associationStatistics.statistics(); }
    public double[] pValues() { return associationStatistics.pValues(); }
    public double[] degreesOfFreedom() {
        return associationStatistics.degreesOfFreedom();
    }
    public double restrictedLogLikelihood() { return restrictedLogLikelihood; }
    /** Maximized profile ML or restricted likelihood, as identified by the options. */
    public double logLikelihood() { return restrictedLogLikelihood; }
    public VarianceEstimation varianceEstimation() { return varianceEstimation; }
    public int observations() { return observations; }
    public int fixedEffectCount() { return fixedEffectCount; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public BackendProvenance backend() { return backend; }

    /** Detects a variance component effectively on the zero boundary. */
    public boolean isSingular(double relativeTolerance) {
        if (!(relativeTolerance > 0.0) || !Double.isFinite(relativeTolerance)) {
            throw new IllegalArgumentException(
                "relativeTolerance must be finite and positive");
        }
        double maximum = 0.0;
        for (double value : varianceComponents) maximum = Math.max(maximum, value);
        if (maximum == 0.0) return true;
        for (double value : varianceComponents) {
            if (value <= maximum * relativeTolerance) return true;
        }
        return false;
    }

    /** Largest absolute fitted log-variance score. */
    public double maximumAbsoluteScore() {
        double maximum = 0.0;
        for (double value : logVarianceScore) maximum = Math.max(maximum, Math.abs(value));
        return maximum;
    }

    /**
     * Tests a fixed-effect contrast using the selected inference covariance.
     * Multi-row tests use the smallest finite coefficient DF involved, a
     * conservative extension of coefficient-level Satterthwaite/KR inference.
     */
    public ContrastTestResult testContrast(double[][] contrast) {
        double denominator = Double.POSITIVE_INFINITY;
        for (double value : degreesOfFreedom()) {
            if (Double.isFinite(value)) denominator = Math.min(denominator, value);
        }
        if (!Double.isFinite(denominator)) {
            return LinearHypothesis.chiSquareTest(
                fixedEffects, fixedEffectInferenceCovariance, contrast);
        }
        return LinearHypothesis.fTest(fixedEffects,
            fixedEffectInferenceCovariance, contrast, denominator);
    }
}
