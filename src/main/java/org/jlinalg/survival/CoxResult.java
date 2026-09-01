/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.List;
import jdistlib.Normal;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Fixed-effect Cox estimates, hazard ratios, inference, and baseline hazard. */
public final class CoxResult {
    private final double[] beta;
    private final double[] covariance;
    private final AssociationStatistics statistics;
    private final double[] confidenceLower;
    private final double[] confidenceUpper;
    private final List<BaselineHazardPoint> baselineHazard;
    private final double logPartialLikelihood;
    private final CoxTies ties;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final BackendProvenance backend;

    CoxResult(
            double[] beta, double[] covariance,
            List<BaselineHazardPoint> baselineHazard,
            double logPartialLikelihood, CoxOptions options,
            int iterations, boolean converged, String convergenceMessage,
            BackendProvenance backend) {
        this.beta = beta.clone();
        this.covariance = covariance.clone();
        double[] standardErrors = new double[beta.length];
        for (int column = 0; column < beta.length; column++)
            standardErrors[column] = Math.sqrt(Math.max(0,
                covariance[column * beta.length + column]));
        statistics = AssociationStatistics.normal(beta, standardErrors);
        double critical = Normal.quantile(
            0.5 + options.confidenceLevel() / 2,
            0, 1, true, false);
        confidenceLower = new double[beta.length];
        confidenceUpper = new double[beta.length];
        for (int column = 0; column < beta.length; column++) {
            confidenceLower[column] = Math.exp(
                beta[column] - critical * standardErrors[column]);
            confidenceUpper[column] = Math.exp(
                beta[column] + critical * standardErrors[column]);
        }
        this.baselineHazard = List.copyOf(baselineHazard);
        this.logPartialLikelihood = logPartialLikelihood;
        this.ties = options.ties();
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
        this.backend = backend;
    }

    public double[] beta() { return beta.clone(); }
    public double[] hazardRatios() {
        return java.util.Arrays.stream(beta).map(Math::exp).toArray();
    }
    public double[] standardErrors() { return statistics.standardErrors(); }
    public double[] zStatistics() { return statistics.statistics(); }
    public double[] pValues() { return statistics.pValues(); }
    public double[] negativeLog10PValues() {
        return statistics.negativeLog10PValues();
    }
    public double[] hazardRatioConfidenceLower() {
        return confidenceLower.clone();
    }
    public double[] hazardRatioConfidenceUpper() {
        return confidenceUpper.clone();
    }
    public double[] covariance() { return covariance.clone(); }
    public AssociationStatistics associationStatistics() { return statistics; }
    public List<BaselineHazardPoint> baselineHazard() { return baselineHazard; }
    public double logPartialLikelihood() { return logPartialLikelihood; }
    public CoxTies ties() { return ties; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public BackendProvenance backend() { return backend; }
}
