/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.ContrastTestResult;
import org.jlinalg.inference.LinearHypothesis;

/** Immutable individual-level two-stage least-squares fit and diagnostics. */
public final class InstrumentalVariableResult {
    private final double[] coefficients;
    private final double[] covariance;
    private final double[] standardErrors;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final double[] fittedValues;
    private final double[] residuals;
    private final double[] instrumentedEndogenous;
    private final int observations;
    private final int exogenousParameters;
    private final int endogenousParameters;
    private final int excludedInstruments;
    private final int structuralRank;
    private final int instrumentRank;
    private final int residualDegreesOfFreedom;
    private final double inferenceDegreesOfFreedom;
    private final int clusters;
    private final double residualSumOfSquares;
    private final double residualVariance;
    private final double projectedDesignConditionNumber;
    private final List<InstrumentStrengthDiagnostic> strengthDiagnostics;
    private final AssociationStatistics associationStatistics;
    private final InstrumentalVariableCovariance covarianceEstimator;
    private final BackendProvenance backend;

    InstrumentalVariableResult(
            double[] coefficients,
            double[] covariance,
            double[] standardErrors,
            double[] confidenceLower,
            double[] confidenceUpper,
            double[] fittedValues,
            double[] residuals,
            double[] instrumentedEndogenous,
            int observations,
            int exogenousParameters,
            int endogenousParameters,
            int excludedInstruments,
            int structuralRank,
            int instrumentRank,
            int residualDegreesOfFreedom,
            double inferenceDegreesOfFreedom,
            int clusters,
            double residualSumOfSquares,
            double residualVariance,
            double projectedDesignConditionNumber,
            List<InstrumentStrengthDiagnostic> strengthDiagnostics,
            AssociationStatistics associationStatistics,
            InstrumentalVariableCovariance covarianceEstimator,
            BackendProvenance backend) {
        this.coefficients = coefficients.clone();
        this.covariance = covariance.clone();
        this.standardErrors = standardErrors.clone();
        this.confidenceLower = confidenceLower.clone();
        this.confidenceUpper = confidenceUpper.clone();
        this.fittedValues = fittedValues.clone();
        this.residuals = residuals.clone();
        this.instrumentedEndogenous = instrumentedEndogenous.clone();
        this.observations = observations;
        this.exogenousParameters = exogenousParameters;
        this.endogenousParameters = endogenousParameters;
        this.excludedInstruments = excludedInstruments;
        this.structuralRank = structuralRank;
        this.instrumentRank = instrumentRank;
        this.residualDegreesOfFreedom = residualDegreesOfFreedom;
        this.inferenceDegreesOfFreedom = inferenceDegreesOfFreedom;
        this.clusters = clusters;
        this.residualSumOfSquares = residualSumOfSquares;
        this.residualVariance = residualVariance;
        this.projectedDesignConditionNumber = projectedDesignConditionNumber;
        this.strengthDiagnostics = List.copyOf(strengthDiagnostics);
        this.associationStatistics = Objects.requireNonNull(
            associationStatistics, "associationStatistics");
        this.covarianceEstimator = Objects.requireNonNull(
            covarianceEstimator, "covarianceEstimator");
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    /** Coefficients ordered as all exogenous columns followed by endogenous columns. */
    public double[] coefficients() { return coefficients.clone(); }

    /** Alias for {@link #coefficients()} used by association APIs. */
    public double[] beta() { return coefficients(); }

    /** Selected row-major coefficient covariance. */
    public double[] covariance() { return covariance.clone(); }

    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] confidenceLower() { return confidenceLower.clone(); }
    public double[] confidenceUpper() { return confidenceUpper.clone(); }

    /** Structural fitted values, using observed rather than instrumented regressors. */
    public double[] fittedValues() { return fittedValues.clone(); }

    /** Structural residuals {@code y - X beta}; these drive every IV covariance. */
    public double[] residuals() { return residuals.clone(); }

    /** First-stage fitted endogenous regressors, returned as observations by columns. */
    public double[][] instrumentedEndogenous() {
        double[][] result = new double[observations][endogenousParameters];
        for (int row = 0; row < observations; row++) {
            System.arraycopy(instrumentedEndogenous,
                row * endogenousParameters, result[row], 0,
                endogenousParameters);
        }
        return result;
    }

    public int observations() { return observations; }
    public int exogenousParameters() { return exogenousParameters; }
    public int endogenousParameters() { return endogenousParameters; }
    public int excludedInstruments() { return excludedInstruments; }
    public int parameters() { return exogenousParameters + endogenousParameters; }
    public int totalInstruments() { return exogenousParameters + excludedInstruments; }
    public int overidentificationDegreesOfFreedom() {
        return excludedInstruments - endogenousParameters;
    }
    public int structuralRank() { return structuralRank; }
    public int instrumentRank() { return instrumentRank; }
    public int residualDegreesOfFreedom() { return residualDegreesOfFreedom; }

    /** Denominator DF used for coefficient inference; infinity denotes normal inference. */
    public double inferenceDegreesOfFreedom() { return inferenceDegreesOfFreedom; }

    /** Independent cluster count, or zero for non-clustered inference. */
    public int clusters() { return clusters; }
    public double residualSumOfSquares() { return residualSumOfSquares; }
    public double residualVariance() { return residualVariance; }

    /**
     * Condition number after scaling projected columns by their observed
     * structural-design norms; this is invariant to regressor unit changes.
     */
    public double projectedDesignConditionNumber() {
        return projectedDesignConditionNumber;
    }

    public List<InstrumentStrengthDiagnostic> strengthDiagnostics() {
        return strengthDiagnostics;
    }

    /** Minimum covariance-matched first-stage effective F statistic. */
    public double minimumFirstStageFStatistic() {
        return strengthDiagnostics.stream()
            .mapToDouble(InstrumentStrengthDiagnostic::effectiveFStatistic)
            .min().orElse(Double.NaN);
    }

    public AssociationStatistics associationStatistics() {
        return associationStatistics;
    }

    public InstrumentalVariableCovariance covarianceEstimator() {
        return covarianceEstimator;
    }

    public BackendProvenance backend() { return backend; }

    /** Wald test of one or more coefficient contrasts using the fit's reference law. */
    public ContrastTestResult testContrast(double[][] contrast) {
        return Double.isFinite(inferenceDegreesOfFreedom)
            ? LinearHypothesis.fTest(coefficients, covariance, contrast,
                inferenceDegreesOfFreedom)
            : LinearHypothesis.chiSquareTest(coefficients, covariance, contrast);
    }
}
