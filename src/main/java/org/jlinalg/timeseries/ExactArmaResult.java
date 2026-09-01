/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import org.jlinalg.compute.BackendProvenance;

/** Exact stationary Gaussian ARMA maximum-likelihood estimates. */
public final class ExactArmaResult {
    private final ArimaOrder order;
    private final double[] autoregressive;
    private final double[] movingAverage;
    private final double mean;
    private final double innovationVariance;
    private final double[] coefficientCovariance;
    private final double[] standardErrors;
    private final double logLikelihood;
    private final double aic;
    private final double bic;
    private final int observedValues;
    private final int seriesCount;
    private final int evaluations;
    private final boolean converged;
    private final BackendProvenance backend;

    ExactArmaResult(ArimaOrder order, double[] ar, double[] ma, double mean,
            double innovationVariance, double[] covariance, double[] standardErrors,
            double logLikelihood, double aic, double bic, int observedValues,
            int seriesCount, int evaluations, boolean converged,
            BackendProvenance backend) {
        this.order = order;
        this.autoregressive = ar.clone();
        this.movingAverage = ma.clone();
        this.mean = mean;
        this.innovationVariance = innovationVariance;
        this.coefficientCovariance = covariance.clone();
        this.standardErrors = standardErrors.clone();
        this.logLikelihood = logLikelihood;
        this.aic = aic;
        this.bic = bic;
        this.observedValues = observedValues;
        this.seriesCount = seriesCount;
        this.evaluations = evaluations;
        this.converged = converged;
        this.backend = backend;
    }
    public ArimaOrder order() { return order; }
    public double[] autoregressive() { return autoregressive.clone(); }
    public double[] movingAverage() { return movingAverage.clone(); }
    public double mean() { return mean; }
    public double innovationVariance() { return innovationVariance; }
    /** Order is AR coefficients, MA coefficients, then mean when included. */
    public double[] coefficientCovariance() { return coefficientCovariance.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double logLikelihood() { return logLikelihood; }
    public double aic() { return aic; }
    public double bic() { return bic; }
    public int observedValues() { return observedValues; }
    public int seriesCount() { return seriesCount; }
    public int functionEvaluations() { return evaluations; }
    public boolean converged() { return converged; }
    public BackendProvenance backend() { return backend; }
}
