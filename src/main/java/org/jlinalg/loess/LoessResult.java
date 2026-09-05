/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.loess;

/** Fitted values, residuals, diagnostics, and reusable LOESS predictor. */
public final class LoessResult {
    private final Loess.Prepared prepared;
    private final double[] response;
    private final double[] priorWeights;
    private final double[] robustnessWeights;
    private final double[] fittedValues;
    private final double[] residuals;
    private final double[] leverage;
    private final double effectiveDegreesOfFreedom;
    private final int iterations;

    LoessResult(Loess.Prepared prepared, double[] response,
            double[] priorWeights, double[] robustnessWeights,
            double[] fittedValues, double[] residuals, double[] leverage,
            double effectiveDegreesOfFreedom, int iterations) {
        this.prepared = prepared;
        this.response = response.clone();
        this.priorWeights = priorWeights.clone();
        this.robustnessWeights = robustnessWeights.clone();
        this.fittedValues = fittedValues.clone();
        this.residuals = residuals.clone();
        this.leverage = leverage.clone();
        this.effectiveDegreesOfFreedom = effectiveDegreesOfFreedom;
        this.iterations = iterations;
    }

    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
    public double[] leverage() { return leverage.clone(); }
    public double[] robustnessWeights() { return robustnessWeights.clone(); }
    public double effectiveDegreesOfFreedom() {
        return effectiveDegreesOfFreedom;
    }
    /** R-compatible alias for the trace of the fitted-point smoother matrix. */
    public double traceHat() { return effectiveDegreesOfFreedom; }
    public int iterations() { return iterations; }
    public LoessOptions options() { return prepared.options(); }

    /** Predicts by direct local fitting, including outside the training range. */
    public double[] predict(double[] x) {
        return prepared.predict(response, priorWeights, robustnessWeights, x);
    }
}
