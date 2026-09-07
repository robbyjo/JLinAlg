/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import org.jlinalg.compute.BackendProvenance;

/** Immutable estimates, inference, fitted parameters, and convergence metadata. */
public final class BetaRegressionResult {
    private final double[] meanCoefficients;
    private final double[] precisionCoefficients;
    private final double[] covariance;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] fittedMeans;
    private final double[] fittedPrecisions;
    private final double logLikelihood;
    private final int observations;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final BetaMeanLink meanLink;
    private final BetaPrecisionLink precisionLink;
    private final BackendProvenance backend;

    BetaRegressionResult(
            double[] meanCoefficients,
            double[] precisionCoefficients,
            double[] covariance,
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] fittedMeans,
            double[] fittedPrecisions,
            double logLikelihood,
            int observations,
            int iterations,
            boolean converged,
            String convergenceMessage,
            BetaMeanLink meanLink,
            BetaPrecisionLink precisionLink,
            BackendProvenance backend) {
        this.meanCoefficients = meanCoefficients.clone();
        this.precisionCoefficients = precisionCoefficients.clone();
        this.covariance = covariance.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.fittedMeans = fittedMeans.clone();
        this.fittedPrecisions = fittedPrecisions.clone();
        this.logLikelihood = logLikelihood;
        this.observations = observations;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
        this.meanLink = meanLink;
        this.precisionLink = precisionLink;
        this.backend = backend;
    }

    public double[] meanCoefficients() { return meanCoefficients.clone(); }
    public double[] precisionCoefficients() {
        return precisionCoefficients.clone();
    }
    /** Covariance of mean coefficients followed by precision coefficients. */
    public double[] covariance() { return covariance.clone(); }
    /** Standard errors in the same block order as {@link #covariance()}. */
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double[] fittedPrecisions() { return fittedPrecisions.clone(); }
    public double logLikelihood() { return logLikelihood; }
    public double aic() {
        return 2.0 * (meanCoefficients.length + precisionCoefficients.length)
            - 2.0 * logLikelihood;
    }
    public int observations() { return observations; }
    public int parameters() {
        return meanCoefficients.length + precisionCoefficients.length;
    }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public BetaMeanLink meanLink() { return meanLink; }
    public BetaPrecisionLink precisionLink() { return precisionLink; }
    public BackendProvenance backend() { return backend; }
}
