/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.reml.VarianceEstimation;

/** Gaussian LMM result with unstructured correlated random blocks. */
public final class CorrelatedLinearMixedModelResult {
    private final AssociationStatistics association;
    private final double[] fixedCovariance;
    private final List<CorrelatedRandomEffectEstimates> randomEffects;
    private final double residualVariance;
    private final double[] fitted;
    private final double[] residuals;
    private final double logLikelihood;
    private final VarianceEstimation estimation;
    private final int functionEvaluations;
    private final boolean converged;
    private final BackendProvenance backend;

    CorrelatedLinearMixedModelResult(
            AssociationStatistics association, double[] fixedCovariance,
            List<CorrelatedRandomEffectEstimates> randomEffects,
            double residualVariance, double[] fitted, double[] residuals,
            double logLikelihood, VarianceEstimation estimation,
            int functionEvaluations, boolean converged,
            BackendProvenance backend) {
        this.association = association;
        this.fixedCovariance = fixedCovariance.clone();
        this.randomEffects = List.copyOf(randomEffects);
        this.residualVariance = residualVariance;
        this.fitted = fitted.clone();
        this.residuals = residuals.clone();
        this.logLikelihood = logLikelihood;
        this.estimation = estimation;
        this.functionEvaluations = functionEvaluations;
        this.converged = converged;
        this.backend = backend;
    }

    public AssociationStatistics associationStatistics() { return association; }
    public double[] beta() { return association.beta(); }
    public double[] fixef() { return association.beta(); }
    public double[] standardErrors() { return association.standardErrors(); }
    public double[] tStatistics() { return association.statistics(); }
    public double[] pValues() { return association.pValues(); }
    public double[] fixedEffectCovariance() { return fixedCovariance.clone(); }
    public List<CorrelatedRandomEffectEstimates> randomEffects() {
        return randomEffects;
    }
    public double residualVariance() { return residualVariance; }
    public double[] fittedValues() { return fitted.clone(); }
    public double[] residuals() { return residuals.clone(); }
    public double logLikelihood() { return logLikelihood; }
    public VarianceEstimation varianceEstimation() { return estimation; }
    public int functionEvaluations() { return functionEvaluations; }
    public boolean converged() { return converged; }
    public BackendProvenance backend() { return backend; }
}
