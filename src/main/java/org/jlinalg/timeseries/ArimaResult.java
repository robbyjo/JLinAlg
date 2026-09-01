/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import java.util.Arrays;

/** Immutable conditional Gaussian ARIMA estimates and diagnostics. */
public final class ArimaResult {
    private final ArimaOrder order;
    private final SeasonalArimaOrder seasonalOrder;
    private final double[] autoregressive;
    private final double[] movingAverage;
    private final double[] seasonalAutoregressive;
    private final double[] seasonalMovingAverage;
    private final double location;
    private final boolean drift;
    private final double innovationVariance;
    private final double[] innovations;
    private final double[] differencedSeries;
    private final double[] originalSeries;
    private final double logLikelihood;
    private final double aic;
    private final double aicc;
    private final double bic;
    private final int effectiveObservations;
    private final int functionEvaluations;
    private final boolean converged;
    private final String convergenceMessage;

    ArimaResult(
            ArimaOrder order,
            SeasonalArimaOrder seasonalOrder,
            double[] autoregressive,
            double[] movingAverage,
            double[] seasonalAutoregressive,
            double[] seasonalMovingAverage,
            double location,
            boolean drift,
            double innovationVariance,
            double[] innovations,
            double[] differencedSeries,
            double[] originalSeries,
            double logLikelihood,
            double aic,
            double aicc,
            double bic,
            int effectiveObservations,
            int functionEvaluations,
            boolean converged,
            String convergenceMessage) {
        this.order = order;
        this.seasonalOrder = seasonalOrder;
        this.autoregressive = autoregressive.clone();
        this.movingAverage = movingAverage.clone();
        this.seasonalAutoregressive = seasonalAutoregressive.clone();
        this.seasonalMovingAverage = seasonalMovingAverage.clone();
        this.location = location;
        this.drift = drift;
        this.innovationVariance = innovationVariance;
        this.innovations = innovations.clone();
        this.differencedSeries = differencedSeries.clone();
        this.originalSeries = originalSeries.clone();
        this.logLikelihood = logLikelihood;
        this.aic = aic;
        this.aicc = aicc;
        this.bic = bic;
        this.effectiveObservations = effectiveObservations;
        this.functionEvaluations = functionEvaluations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
    }

    public ArimaOrder order() { return order; }
    public SeasonalArimaOrder seasonalOrder() { return seasonalOrder; }
    public double[] autoregressive() { return autoregressive.clone(); }
    public double[] movingAverage() { return movingAverage.clone(); }
    public double[] seasonalAutoregressive() {
        return seasonalAutoregressive.clone();
    }
    public double[] seasonalMovingAverage() {
        return seasonalMovingAverage.clone();
    }
    /** Mean for undifferenced ARMA, or drift after one total difference. */
    public double location() { return location; }
    public boolean drift() { return drift; }
    public double innovationVariance() { return innovationVariance; }
    public double[] innovations() { return innovations.clone(); }
    public double[] differencedSeries() { return differencedSeries.clone(); }
    public double logLikelihood() { return logLikelihood; }
    public double aic() { return aic; }
    public double aicc() { return aicc; }
    public double bic() { return bic; }
    public int effectiveObservations() { return effectiveObservations; }
    public int functionEvaluations() { return functionEvaluations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }

    public ArimaForecast forecast(int horizon) {
        return forecast(horizon, 0.95);
    }

    public ArimaForecast forecast(int horizon, double confidenceLevel) {
        return ArimaMath.forecast(this, originalSeries, horizon, confidenceLevel);
    }

    /** Residual autocorrelations after the conditional startup interval. */
    public double[] residualAutocorrelation(int maximumLag) {
        return TimeSeriesDiagnostics.autocorrelation(usedInnovations(), maximumLag);
    }

    /** Ljung-Box diagnostic after accounting for fitted AR and MA parameters. */
    public LjungBoxResult ljungBox(int lags) {
        int fitted = order.autoregressive() + order.movingAverage()
            + seasonalOrder.autoregressive() + seasonalOrder.movingAverage();
        return TimeSeriesDiagnostics.ljungBox(usedInnovations(), lags, fitted);
    }

    private double[] usedInnovations() {
        return Arrays.copyOfRange(innovations,
            innovations.length - effectiveObservations, innovations.length);
    }
}
