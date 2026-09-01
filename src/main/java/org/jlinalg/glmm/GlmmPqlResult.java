/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glmm;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Immutable conditional estimates and convergence metadata from PQL. */
public final class GlmmPqlResult {
    private final String family;
    private final List<String> componentNames;
    private final double[] varianceComponents;
    private final double[] fixedEffects;
    private final double[] fixedEffectCovariance;
    private final double[] fixedEffectInferenceCovariance;
    private final double[] fixedEffectStandardErrors;
    private final AssociationStatistics associationStatistics;
    private final double[] randomLinearPredictor;
    private final double[] linearPredictor;
    private final double[] fittedMeans;
    private final double conditionalDeviance;
    private final double finalWorkingRestrictedLogLikelihood;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final BackendProvenance backend;

    GlmmPqlResult(
            String family, List<String> componentNames,
            double[] varianceComponents, double[] fixedEffects,
            double[] fixedEffectCovariance,
            double[] fixedEffectInferenceCovariance,
            double[] fixedEffectStandardErrors,
            AssociationStatistics associationStatistics,
            double[] randomLinearPredictor, double[] linearPredictor,
            double[] fittedMeans, double conditionalDeviance,
            double finalWorkingRestrictedLogLikelihood,
            int iterations, boolean converged, String convergenceMessage,
            BackendProvenance backend) {
        this.family = Objects.requireNonNull(family, "family");
        this.componentNames = List.copyOf(componentNames);
        this.varianceComponents = varianceComponents.clone();
        this.fixedEffects = fixedEffects.clone();
        this.fixedEffectCovariance = fixedEffectCovariance.clone();
        this.fixedEffectInferenceCovariance =
            fixedEffectInferenceCovariance.clone();
        this.fixedEffectStandardErrors = fixedEffectStandardErrors.clone();
        this.associationStatistics = Objects.requireNonNull(
            associationStatistics, "associationStatistics");
        this.randomLinearPredictor = randomLinearPredictor.clone();
        this.linearPredictor = linearPredictor.clone();
        this.fittedMeans = fittedMeans.clone();
        this.conditionalDeviance = conditionalDeviance;
        this.finalWorkingRestrictedLogLikelihood =
            finalWorkingRestrictedLogLikelihood;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = Objects.requireNonNull(
            convergenceMessage, "convergenceMessage");
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    public String family() { return family; }
    public List<String> componentNames() { return componentNames; }
    public double[] varianceComponents() { return varianceComponents.clone(); }
    public double[] fixedEffects() { return fixedEffects.clone(); }
    /** Unadjusted covariance from the final working REML model. */
    public double[] fixedEffectCovariance() { return fixedEffectCovariance.clone(); }
    /** Covariance used for association inference; KR adjusts this matrix. */
    public double[] fixedEffectInferenceCovariance() {
        return fixedEffectInferenceCovariance.clone();
    }
    public double[] fixedEffectStandardErrors() {
        return fixedEffectStandardErrors.clone();
    }
    /** Returns beta, SE, t statistics, denominator DF, and two-sided p-values. */
    public AssociationStatistics associationStatistics() {
        return associationStatistics;
    }
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
    public double[] randomLinearPredictor() { return randomLinearPredictor.clone(); }
    public double[] linearPredictor() { return linearPredictor.clone(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double conditionalDeviance() { return conditionalDeviance; }
    public double finalWorkingRestrictedLogLikelihood() {
        return finalWorkingRestrictedLogLikelihood;
    }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public BackendProvenance backend() { return backend; }
}
