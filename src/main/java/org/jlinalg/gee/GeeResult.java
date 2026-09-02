/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.ContrastTestResult;
import org.jlinalg.inference.LinearHypothesis;

/** Immutable GEE estimates, covariance variants, criteria, and diagnostics. */
public final class GeeResult {
    private final String family;
    private final GeeCorrelation correlation;
    private final GeeAssociation association;
    private final GeeCovariance covarianceType;
    private final GeeMethod method;
    private final double[] coefficients;
    private final double[] covariance;
    private final double[] naiveCovariance;
    private final double[] robustCovariance;
    private final double[] dfAdjustedCovariance;
    private final double[] biasCorrectedCovariance;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final double[] linearPredictor;
    private final double[] fittedMeans;
    private final double[] pearsonResiduals;
    private final double[] associationParameters;
    private final double dispersion;
    private final double[] scaleCoefficients;
    private final GeeCriteria criteria;
    private final int observations;
    private final int clusters;
    private final int minimumClusterSize;
    private final int maximumClusterSize;
    private final int parameters;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final int[] retainedRows;
    private final int originalObservations;
    private final BackendProvenance backend;

    GeeResult(
            String family,
            GeeCorrelation correlation,
            GeeAssociation association,
            GeeCovariance covarianceType,
            GeeMethod method,
            double[] coefficients,
            double[] covariance,
            double[] naiveCovariance,
            double[] robustCovariance,
            double[] dfAdjustedCovariance,
            double[] biasCorrectedCovariance,
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] confidenceLower,
            double[] confidenceUpper,
            double[] linearPredictor,
            double[] fittedMeans,
            double[] pearsonResiduals,
            double[] associationParameters,
            double dispersion,
            double[] scaleCoefficients,
            GeeCriteria criteria,
            int observations,
            int clusters,
            int minimumClusterSize,
            int maximumClusterSize,
            int parameters,
            int iterations,
            boolean converged,
            String convergenceMessage,
            int[] retainedRows,
            int originalObservations,
            BackendProvenance backend) {
        this.family = Objects.requireNonNull(family, "family");
        this.correlation = Objects.requireNonNull(correlation, "correlation");
        this.association = Objects.requireNonNull(association, "association");
        this.covarianceType = Objects.requireNonNull(covarianceType, "covarianceType");
        this.method = Objects.requireNonNull(method, "method");
        this.coefficients = coefficients.clone();
        this.covariance = covariance.clone();
        this.naiveCovariance = naiveCovariance.clone();
        this.robustCovariance = robustCovariance.clone();
        this.dfAdjustedCovariance = dfAdjustedCovariance.clone();
        this.biasCorrectedCovariance = biasCorrectedCovariance.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.confidenceLower = confidenceLower.clone();
        this.confidenceUpper = confidenceUpper.clone();
        this.linearPredictor = linearPredictor.clone();
        this.fittedMeans = fittedMeans.clone();
        this.pearsonResiduals = pearsonResiduals.clone();
        this.associationParameters = associationParameters.clone();
        this.dispersion = dispersion;
        this.scaleCoefficients = scaleCoefficients.clone();
        this.criteria = Objects.requireNonNull(criteria, "criteria");
        this.observations = observations;
        this.clusters = clusters;
        this.minimumClusterSize = minimumClusterSize;
        this.maximumClusterSize = maximumClusterSize;
        this.parameters = parameters;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = Objects.requireNonNull(
            convergenceMessage, "convergenceMessage");
        this.retainedRows = retainedRows.clone();
        this.originalObservations = originalObservations;
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    public String family() { return family; }
    public GeeCorrelation correlation() { return correlation; }
    public GeeAssociation association() { return association; }
    public GeeCovariance covarianceType() { return covarianceType; }
    public GeeMethod method() { return method; }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] beta() { return coefficients(); }
    public double[] covariance() { return covariance.clone(); }
    public double[] naiveCovariance() { return naiveCovariance.clone(); }
    public double[] robustCovariance() { return robustCovariance.clone(); }
    public double[] dfAdjustedCovariance() { return dfAdjustedCovariance.clone(); }
    public double[] biasCorrectedCovariance() {
        return biasCorrectedCovariance.clone();
    }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] confidenceLower() { return confidenceLower.clone(); }
    public double[] confidenceUpper() { return confidenceUpper.clone(); }
    public double[] linearPredictor() { return linearPredictor.clone(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double[] pearsonResiduals() { return pearsonResiduals.clone(); }
    public double[] associationParameters() {
        return associationParameters.clone();
    }
    public double dispersion() { return dispersion; }
    public double[] scaleCoefficients() { return scaleCoefficients.clone(); }
    public GeeCriteria criteria() { return criteria; }
    public int observations() { return observations; }
    public int clusters() { return clusters; }
    public int minimumClusterSize() { return minimumClusterSize; }
    public int maximumClusterSize() { return maximumClusterSize; }
    public int parameters() { return parameters; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public int[] retainedRows() { return retainedRows.clone(); }
    public int originalObservations() { return originalObservations; }
    public int omittedObservations() { return originalObservations - observations; }
    public BackendProvenance backend() { return backend; }

    /** Returns coefficient-level asymptotic Wald z inference. */
    public AssociationStatistics associationStatistics() {
        return AssociationStatistics.normal(coefficients, standardErrors);
    }

    /** Tests one or more linear contrasts using the selected covariance. */
    public ContrastTestResult testContrast(double[][] contrast) {
        return LinearHypothesis.chiSquareTest(coefficients, covariance, contrast);
    }
}
