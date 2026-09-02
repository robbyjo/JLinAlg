/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Objects;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.ContrastTestResult;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.inference.LinearHypothesis;
import org.jlinalg.glm.GlmFamily;

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
    private final double[] kauermannCarrollCovariance;
    private final double[] fayGraubardCovariance;
    private final double[] jackknifeCovariance;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final double[] linearPredictor;
    private final double[] fittedMeans;
    private final double[] pearsonResiduals;
    private final double[] responseResiduals;
    private final double[] devianceResiduals;
    private final double[] workingResiduals;
    private final double[] standardizedResiduals;
    private final double[] associationParameters;
    private final double dispersion;
    private final double[] scaleCoefficients;
    private final GeeCriteria criteria;
    private final GeeInference inference;
    private final double degreesOfFreedom;
    private final GeeDiagnostics diagnostics;
    private final GeeConvergenceDiagnostics convergenceDiagnostics;
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
            double[] kauermannCarrollCovariance,
            double[] fayGraubardCovariance,
            double[] jackknifeCovariance,
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] confidenceLower,
            double[] confidenceUpper,
            double[] linearPredictor,
            double[] fittedMeans,
            double[] pearsonResiduals,
            double[] responseResiduals,
            double[] devianceResiduals,
            double[] workingResiduals,
            double[] standardizedResiduals,
            double[] associationParameters,
            double dispersion,
            double[] scaleCoefficients,
            GeeCriteria criteria,
            GeeInference inference,
            double degreesOfFreedom,
            GeeDiagnostics diagnostics,
            GeeConvergenceDiagnostics convergenceDiagnostics,
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
        this.kauermannCarrollCovariance = kauermannCarrollCovariance.clone();
        this.fayGraubardCovariance = fayGraubardCovariance.clone();
        this.jackknifeCovariance = jackknifeCovariance.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.confidenceLower = confidenceLower.clone();
        this.confidenceUpper = confidenceUpper.clone();
        this.linearPredictor = linearPredictor.clone();
        this.fittedMeans = fittedMeans.clone();
        this.pearsonResiduals = pearsonResiduals.clone();
        this.responseResiduals = responseResiduals.clone();
        this.devianceResiduals = devianceResiduals.clone();
        this.workingResiduals = workingResiduals.clone();
        this.standardizedResiduals = standardizedResiduals.clone();
        this.associationParameters = associationParameters.clone();
        this.dispersion = dispersion;
        this.scaleCoefficients = scaleCoefficients.clone();
        this.criteria = Objects.requireNonNull(criteria, "criteria");
        this.inference = Objects.requireNonNull(inference, "inference");
        this.degreesOfFreedom = degreesOfFreedom;
        this.diagnostics = Objects.requireNonNull(diagnostics, "diagnostics");
        this.convergenceDiagnostics = Objects.requireNonNull(
            convergenceDiagnostics, "convergenceDiagnostics");
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
    public double[] kauermannCarrollCovariance() {
        return kauermannCarrollCovariance.clone();
    }
    public double[] fayGraubardCovariance() {
        return fayGraubardCovariance.clone();
    }
    public double[] jackknifeCovariance() { return jackknifeCovariance.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] confidenceLower() { return confidenceLower.clone(); }
    public double[] confidenceUpper() { return confidenceUpper.clone(); }
    public double[] linearPredictor() { return linearPredictor.clone(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double[] pearsonResiduals() { return pearsonResiduals.clone(); }
    public double[] responseResiduals() { return responseResiduals.clone(); }
    public double[] devianceResiduals() { return devianceResiduals.clone(); }
    public double[] workingResiduals() { return workingResiduals.clone(); }
    public double[] standardizedResiduals() {
        return standardizedResiduals.clone();
    }
    public double[] residuals(GeeResidualType type) {
        Objects.requireNonNull(type, "type");
        return switch (type) {
            case RESPONSE -> responseResiduals();
            case PEARSON -> pearsonResiduals();
            case DEVIANCE -> devianceResiduals();
            case WORKING -> workingResiduals();
            case STANDARDIZED -> standardizedResiduals();
        };
    }
    public double[] associationParameters() {
        return associationParameters.clone();
    }
    public double dispersion() { return dispersion; }
    public double[] scaleCoefficients() { return scaleCoefficients.clone(); }
    public GeeCriteria criteria() { return criteria; }
    public GeeInference inference() { return inference; }
    public double degreesOfFreedom() { return degreesOfFreedom; }
    public GeeDiagnostics diagnostics() { return diagnostics; }
    public GeeConvergenceDiagnostics convergenceDiagnostics() {
        return convergenceDiagnostics;
    }
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
        return inference == GeeInference.CLUSTER_T
            ? AssociationStatistics.studentT(coefficients, standardErrors,
                degreesOfFreedom, DegreesOfFreedomMethod.CLUSTER)
            : AssociationStatistics.normal(coefficients, standardErrors);
    }

    /** Tests one or more linear contrasts using the selected covariance. */
    public ContrastTestResult testContrast(double[][] contrast) {
        return inference == GeeInference.CLUSTER_T
            ? LinearHypothesis.fTest(coefficients, covariance, contrast,
                degreesOfFreedom)
            : LinearHypothesis.chiSquareTest(coefficients, covariance, contrast);
    }

    /** Returns one immutable row per coefficient for reporting or export. */
    public GeeCoefficient[] coefficientTable() {
        GeeCoefficient[] result = new GeeCoefficient[coefficients.length];
        for (int index = 0; index < result.length; index++) {
            result[index] = new GeeCoefficient(index, coefficients[index],
                standardErrors[index], statistics[index], pValues[index],
                confidenceLower[index], confidenceUpper[index], degreesOfFreedom);
        }
        return result;
    }

    /** Serializes the principal estimates and reproducibility metadata as JSON. */
    public String toJson() {
        StringBuilder json = new StringBuilder(512);
        json.append('{');
        appendString(json, "family", family).append(',');
        appendString(json, "correlation", correlation.name()).append(',');
        appendString(json, "association", association.name()).append(',');
        appendString(json, "covarianceType", covarianceType.name()).append(',');
        appendString(json, "method", method.name()).append(',');
        appendString(json, "inference", inference.name()).append(',');
        appendArray(json, "coefficients", coefficients).append(',');
        appendArray(json, "covariance", covariance).append(',');
        appendArray(json, "standardErrors", standardErrors).append(',');
        appendArray(json, "statistics", statistics).append(',');
        appendArray(json, "pValues", pValues).append(',');
        appendArray(json, "confidenceLower", confidenceLower).append(',');
        appendArray(json, "confidenceUpper", confidenceUpper).append(',');
        appendArray(json, "associationParameters", associationParameters)
            .append(',');
        appendNumber(json, "dispersion", dispersion).append(',');
        appendArray(json, "scaleCoefficients", scaleCoefficients).append(',');
        appendNumber(json, "quasiLikelihood", criteria.quasiLikelihood())
            .append(',');
        appendNumber(json, "qic", criteria.qic()).append(',');
        appendNumber(json, "qicu", criteria.qicu()).append(',');
        appendNumber(json, "cic", criteria.cic()).append(',');
        appendNumber(json, "qicc", criteria.qicc()).append(',');
        appendNumber(json, "degreesOfFreedom", degreesOfFreedom).append(',');
        appendInteger(json, "observations", observations).append(',');
        appendInteger(json, "clusters", clusters).append(',');
        appendInteger(json, "parameters", parameters).append(',');
        appendInteger(json, "iterations", iterations).append(',');
        json.append("\"converged\":").append(converged).append(',');
        appendString(json, "convergenceMessage", convergenceMessage).append(',');
        appendString(json, "requestedBackend", backend.requested().name())
            .append(',');
        appendString(json, "selectedBackend", backend.selectedBackend())
            .append(',');
        appendString(json, "deviceDescription", backend.deviceDescription());
        return json.append('}').toString();
    }

    /** Predicts marginal means for new design rows using the selected covariance. */
    public GeePrediction[] predict(double[][] design, GlmFamily family) {
        return ClinicalGee.predict(this, design, family);
    }

    /** Predicts marginal means for new design rows and supplied offsets. */
    public GeePrediction[] predict(
            double[][] design,
            GlmFamily family,
            double[] offset,
            double confidenceLevel) {
        return ClinicalGee.predict(this, design, family, offset, confidenceLevel);
    }

    private static StringBuilder appendString(
            StringBuilder json, String name, String value) {
        appendName(json, name).append('"');
        for (int index = 0; index < value.length(); index++) {
            char character = value.charAt(index);
            switch (character) {
                case '"' -> json.append("\\\"");
                case '\\' -> json.append("\\\\");
                case '\b' -> json.append("\\b");
                case '\f' -> json.append("\\f");
                case '\n' -> json.append("\\n");
                case '\r' -> json.append("\\r");
                case '\t' -> json.append("\\t");
                default -> {
                    if (character < 0x20) {
                        json.append(String.format("\\u%04x", (int) character));
                    } else {
                        json.append(character);
                    }
                }
            }
        }
        return json.append('"');
    }

    private static StringBuilder appendArray(
            StringBuilder json, String name, double[] values) {
        appendName(json, name).append('[');
        for (int index = 0; index < values.length; index++) {
            if (index > 0) json.append(',');
            appendValue(json, values[index]);
        }
        return json.append(']');
    }

    private static StringBuilder appendNumber(
            StringBuilder json, String name, double value) {
        appendName(json, name);
        return appendValue(json, value);
    }

    private static StringBuilder appendInteger(
            StringBuilder json, String name, int value) {
        return appendName(json, name).append(value);
    }

    private static StringBuilder appendName(StringBuilder json, String name) {
        return json.append('"').append(name).append("\":");
    }

    private static StringBuilder appendValue(StringBuilder json, double value) {
        return Double.isFinite(value) ? json.append(value) : json.append("null");
    }
}
