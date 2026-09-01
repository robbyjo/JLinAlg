/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.timeseries;

import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.mixed.LinearMixedModelResult;

/** Profile-REML LMM estimates with ARIMA-correlated errors. */
public final class ArimaErrorLmmResult {
    private final LinearMixedModelResult mixedModel;
    private final ArimaOrder order;
    private final SeasonalArimaOrder seasonalOrder;
    private final double[] autoregressive;
    private final double[] movingAverage;
    private final double[] seasonalAutoregressive;
    private final double[] seasonalMovingAverage;
    private final int differencingLoss;
    private final int functionEvaluations;
    private final boolean converged;
    private final String convergenceMessage;

    ArimaErrorLmmResult(
            LinearMixedModelResult mixedModel,
            ArimaOrder order,
            SeasonalArimaOrder seasonalOrder,
            ArimaMath.Coefficients coefficients,
            int differencingLoss,
            int functionEvaluations,
            boolean converged,
            String convergenceMessage) {
        this.mixedModel = mixedModel;
        this.order = order;
        this.seasonalOrder = seasonalOrder;
        this.autoregressive = coefficients.ar().clone();
        this.movingAverage = coefficients.ma().clone();
        this.seasonalAutoregressive = coefficients.seasonalAr().clone();
        this.seasonalMovingAverage = coefficients.seasonalMa().clone();
        this.differencingLoss = differencingLoss;
        this.functionEvaluations = functionEvaluations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
    }

    public LinearMixedModelResult mixedModel() { return mixedModel; }
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
    public int differencingLoss() { return differencingLoss; }
    public boolean differenced() { return differencingLoss > 0; }
    public int functionEvaluations() { return functionEvaluations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public AssociationStatistics associationStatistics() {
        return mixedModel.associationStatistics();
    }
    public double[] beta() { return mixedModel.beta(); }
    public double[] standardErrors() { return mixedModel.standardErrors(); }
    public double[] tStatistics() { return mixedModel.tStatistics(); }
    public double[] pValues() { return mixedModel.pValues(); }
}
