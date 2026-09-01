/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.penalized;

import org.jlinalg.inference.AssociationStatistics;

/** Ridge fit with effective-DF model-based coefficient inference. */
public final class RidgeRegressionResult {
    private final PenalizedRegressionResult fit;
    private final AssociationStatistics associationStatistics;
    private final double[] coefficientCovariance;
    private final double effectiveModelDegreesOfFreedom;
    private final double residualDegreesOfFreedom;
    private final double residualVariance;

    RidgeRegressionResult(
            PenalizedRegressionResult fit,
            AssociationStatistics associationStatistics,
            double[] coefficientCovariance,
            double effectiveModelDegreesOfFreedom,
            double residualDegreesOfFreedom,
            double residualVariance) {
        this.fit = fit;
        this.associationStatistics = associationStatistics;
        this.coefficientCovariance = coefficientCovariance.clone();
        this.effectiveModelDegreesOfFreedom = effectiveModelDegreesOfFreedom;
        this.residualDegreesOfFreedom = residualDegreesOfFreedom;
        this.residualVariance = residualVariance;
    }

    public PenalizedRegressionResult fit() { return fit; }
    public AssociationStatistics associationStatistics() {
        return associationStatistics;
    }
    public double intercept() { return fit.intercept(); }
    public double[] beta() { return fit.beta(); }
    public double[] standardErrors() {
        return associationStatistics.standardErrors();
    }
    public double[] tStatistics() { return associationStatistics.statistics(); }
    public double[] pValues() { return associationStatistics.pValues(); }
    public double[] coefficientCovariance() {
        return coefficientCovariance.clone();
    }
    public double effectiveModelDegreesOfFreedom() {
        return effectiveModelDegreesOfFreedom;
    }
    public double residualDegreesOfFreedom() {
        return residualDegreesOfFreedom;
    }
    public double residualVariance() { return residualVariance; }
}
