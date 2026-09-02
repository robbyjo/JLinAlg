/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.ArrayList;
import java.util.List;

/** Penalized Gaussian fit and selected multi-penalty smoothing parameters. */
public final class GaussianSmoothSelectionResult {
    private final PenalizedPredictor predictor;
    private final List<double[]> smoothingParameters;
    private final double[] coefficients;
    private final double[] fittedValues;
    private final double[] residuals;
    private final double covarianceScale;
    private final double[] covariance;
    private final double effectiveDegreesOfFreedom;
    private final double criterionValue;
    private final int evaluations;

    GaussianSmoothSelectionResult(
            PenalizedPredictor predictor,
            List<double[]> smoothingParameters,
            double[] coefficients,
            double[] fittedValues,
            double[] residuals,
            double covarianceScale,
            double[] covariance,
            double effectiveDegreesOfFreedom,
            double criterionValue,
            int evaluations) {
        this.predictor = predictor;
        List<double[]> copied = new ArrayList<>(smoothingParameters.size());
        for (double[] values : smoothingParameters) copied.add(values.clone());
        this.smoothingParameters = List.copyOf(copied);
        this.coefficients = coefficients.clone();
        this.fittedValues = fittedValues.clone();
        this.residuals = residuals.clone();
        this.covarianceScale = covarianceScale;
        this.covariance = covariance.clone();
        this.effectiveDegreesOfFreedom = effectiveDegreesOfFreedom;
        this.criterionValue = criterionValue;
        this.evaluations = evaluations;
    }

    public PenalizedPredictor predictor() { return predictor; }
    public List<double[]> smoothingParameters() {
        List<double[]> result = new ArrayList<>(smoothingParameters.size());
        for (double[] values : smoothingParameters) result.add(values.clone());
        return List.copyOf(result);
    }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double[] residuals() { return residuals.clone(); }
    public double residualVariance() { return covarianceScale; }
    public double[] covariance() { return covariance.clone(); }
    public double effectiveDegreesOfFreedom() { return effectiveDegreesOfFreedom; }
    public double criterionValue() { return criterionValue; }
    public int evaluations() { return evaluations; }
}
