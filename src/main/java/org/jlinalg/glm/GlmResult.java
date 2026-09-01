/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.ContrastTestResult;
import org.jlinalg.inference.LinearHypothesis;

/** Immutable GLM estimates, diagnostics, inference, and convergence metadata. */
public final class GlmResult {
    private final String family;
    private final double[] coefficients;
    private final double[] covariance;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final double[] linearPredictor;
    private final double[] fittedMeans;
    private final double[] devianceResiduals;
    private final double[] pearsonResiduals;
    private final double deviance;
    private final double dispersion;
    private final double logLikelihood;
    private final double aic;
    private final int observations;
    private final int parameters;
    private final int rank;
    private final int residualDegreesOfFreedom;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final int[] retainedRows;
    private final int originalObservations;
    private final BackendProvenance backend;

    GlmResult(
            String family, double[] coefficients, double[] covariance,
            double[] standardErrors, double[] statistics, double[] pValues,
            double[] confidenceLower, double[] confidenceUpper,
            double[] linearPredictor, double[] fittedMeans,
            double[] devianceResiduals, double[] pearsonResiduals,
            double deviance, double dispersion, double logLikelihood, double aic,
            int observations, int parameters, int rank,
            int residualDegreesOfFreedom, int iterations,
            boolean converged, String convergenceMessage,
            int[] retainedRows, int originalObservations,
            BackendProvenance backend) {
        this.family = Objects.requireNonNull(family, "family");
        this.coefficients = coefficients.clone();
        this.covariance = covariance.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.confidenceLower = confidenceLower.clone();
        this.confidenceUpper = confidenceUpper.clone();
        this.linearPredictor = linearPredictor.clone();
        this.fittedMeans = fittedMeans.clone();
        this.devianceResiduals = devianceResiduals.clone();
        this.pearsonResiduals = pearsonResiduals.clone();
        this.deviance = deviance;
        this.dispersion = dispersion;
        this.logLikelihood = logLikelihood;
        this.aic = aic;
        this.observations = observations;
        this.parameters = parameters;
        this.rank = rank;
        this.residualDegreesOfFreedom = residualDegreesOfFreedom;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = Objects.requireNonNull(
            convergenceMessage, "convergenceMessage");
        this.retainedRows = retainedRows.clone();
        this.originalObservations = originalObservations;
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    public String family() { return family; }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] covariance() { return covariance.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] confidenceLower() { return confidenceLower.clone(); }
    public double[] confidenceUpper() { return confidenceUpper.clone(); }
    public double[] linearPredictor() { return linearPredictor.clone(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double[] devianceResiduals() { return devianceResiduals.clone(); }
    public double[] pearsonResiduals() { return pearsonResiduals.clone(); }
    public double deviance() { return deviance; }
    public double dispersion() { return dispersion; }
    public double logLikelihood() { return logLikelihood; }
    public double aic() { return aic; }
    public int observations() { return observations; }
    public int parameters() { return parameters; }
    public int rank() { return rank; }
    public int residualDegreesOfFreedom() { return residualDegreesOfFreedom; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    /** Original zero-based row indices used in the fit. */
    public int[] retainedRows() { return retainedRows.clone(); }
    public int originalObservations() { return originalObservations; }
    public int omittedObservations() { return originalObservations - observations; }
    public BackendProvenance backend() { return backend; }

    /** Returns coefficient-level asymptotic Wald z inference. */
    public AssociationStatistics associationStatistics() {
        return AssociationStatistics.normal(coefficients, standardErrors);
    }

    /** Alias for {@link #coefficients()} for consistent association APIs. */
    public double[] beta() { return coefficients(); }

    /** Tests one or more linear contrasts using an asymptotic Wald chi-square. */
    public ContrastTestResult testContrast(double[][] contrast) {
        return LinearHypothesis.chiSquareTest(coefficients, covariance, contrast);
    }
}
