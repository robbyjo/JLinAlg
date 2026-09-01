/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.ols;

import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.inference.ContrastTestResult;
import org.jlinalg.inference.LinearHypothesis;

/** Immutable OLS estimates, uncertainty, residuals, and numerical metadata. */
public final class OlsResult {
    private final double[] coefficients;
    private final double[] fittedValues;
    private final double[] residuals;
    private final double[] covariance;
    private final double[] standardErrors;
    private final double[] tStatistics;
    private final double[] pValues;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final int observations;
    private final int parameters;
    private final int rank;
    private final int residualDegreesOfFreedom;
    private final double residualSumOfSquares;
    private final double residualVariance;
    private final double logLikelihood;
    private final boolean minimumNorm;
    private final double singularValueTolerance;
    private final int[] retainedRows;
    private final int originalObservations;
    private final BackendProvenance backend;

    OlsResult(
            double[] coefficients,
            double[] fittedValues,
            double[] residuals,
            double[] covariance,
            double[] standardErrors,
            double[] tStatistics,
            double[] pValues,
            double[] confidenceLower,
            double[] confidenceUpper,
            int observations,
            int parameters,
            int rank,
            int residualDegreesOfFreedom,
            double residualSumOfSquares,
            double residualVariance,
            double logLikelihood,
            boolean minimumNorm,
            double singularValueTolerance,
            int[] retainedRows,
            int originalObservations,
            BackendProvenance backend) {
        this.coefficients = coefficients.clone();
        this.fittedValues = fittedValues.clone();
        this.residuals = residuals.clone();
        this.covariance = covariance.clone();
        this.standardErrors = standardErrors.clone();
        this.tStatistics = tStatistics.clone();
        this.pValues = pValues.clone();
        this.confidenceLower = confidenceLower.clone();
        this.confidenceUpper = confidenceUpper.clone();
        this.observations = observations;
        this.parameters = parameters;
        this.rank = rank;
        this.residualDegreesOfFreedom = residualDegreesOfFreedom;
        this.residualSumOfSquares = residualSumOfSquares;
        this.residualVariance = residualVariance;
        this.logLikelihood = logLikelihood;
        this.minimumNorm = minimumNorm;
        this.singularValueTolerance = singularValueTolerance;
        this.retainedRows = retainedRows.clone();
        this.originalObservations = originalObservations;
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    public double[] coefficients() { return coefficients.clone(); }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
    public double[] covariance() { return covariance.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] tStatistics() { return tStatistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] confidenceLower() { return confidenceLower.clone(); }
    public double[] confidenceUpper() { return confidenceUpper.clone(); }
    public int observations() { return observations; }
    public int parameters() { return parameters; }
    public int rank() { return rank; }
    public int residualDegreesOfFreedom() { return residualDegreesOfFreedom; }
    public double residualSumOfSquares() { return residualSumOfSquares; }
    public double residualVariance() { return residualVariance; }
    public double logLikelihood() { return logLikelihood; }
    public boolean rankDeficient() { return rank < parameters; }
    public boolean minimumNorm() { return minimumNorm; }
    public double singularValueTolerance() { return singularValueTolerance; }
    /** Original zero-based row indices used in the fit. */
    public int[] retainedRows() { return retainedRows.clone(); }
    public int originalObservations() { return originalObservations; }
    public int omittedObservations() { return originalObservations - observations; }
    public BackendProvenance backend() { return backend; }

    /** Returns the coefficient-level OLS association table. */
    public AssociationStatistics associationStatistics() {
        return AssociationStatistics.studentT(coefficients, standardErrors,
            residualDegreesOfFreedom,
            DegreesOfFreedomMethod.RESIDUAL);
    }

    /** Alias for {@link #coefficients()} for consistent association APIs. */
    public double[] beta() { return coefficients(); }

    /** Tests one or more estimable linear contrasts using the residual F law. */
    public ContrastTestResult testContrast(double[][] contrast) {
        return LinearHypothesis.fTest(coefficients, covariance,
            contrast, residualDegreesOfFreedom);
    }
}
