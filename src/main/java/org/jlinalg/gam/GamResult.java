/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import org.jlinalg.mixed.LinearMixedModelResult;

/** Exact Gaussian REML GAM result with per-smooth estimates. */
public final class GamResult {
    private final LinearMixedModelResult mixedModel;
    private final int parametricColumns;
    private final List<SmoothTermEstimate> smoothTerms;
    private final double totalEffectiveDegreesOfFreedom;

    GamResult(
            LinearMixedModelResult mixedModel,
            int parametricColumns,
            List<SmoothTermEstimate> smoothTerms,
            double totalEffectiveDegreesOfFreedom) {
        this.mixedModel = mixedModel;
        this.parametricColumns = parametricColumns;
        this.smoothTerms = List.copyOf(smoothTerms);
        this.totalEffectiveDegreesOfFreedom = totalEffectiveDegreesOfFreedom;
    }

    /** Coefficients for the caller-supplied parametric design columns. */
    public double[] parametricCoefficients() {
        return prefix(mixedModel.beta(), parametricColumns);
    }

    /** Standard errors for the caller-supplied parametric design columns. */
    public double[] parametricStandardErrors() {
        return prefix(mixedModel.standardErrors(), parametricColumns);
    }

    /** t statistics for the caller-supplied parametric design columns. */
    public double[] parametricTStatistics() {
        return prefix(mixedModel.tStatistics(), parametricColumns);
    }

    /** p values for the caller-supplied parametric design columns. */
    public double[] parametricPValues() {
        return prefix(mixedModel.pValues(), parametricColumns);
    }

    public List<SmoothTermEstimate> smoothTerms() { return smoothTerms; }
    public double totalEffectiveDegreesOfFreedom() {
        return totalEffectiveDegreesOfFreedom;
    }
    public double residualDegreesOfFreedom() {
        return mixedModel.reml().observations() - totalEffectiveDegreesOfFreedom;
    }
    public double[] fittedValues() { return mixedModel.fittedValues(); }
    public double[] residuals() { return mixedModel.residuals(); }
    public LinearMixedModelResult mixedModel() { return mixedModel; }

    /**
     * Predicts from a row-major parametric design and covariate vectors in the
     * same order as {@link #smoothTerms()}.
     */
    public double[] predict(
            double[] parametricDesign,
            int rows,
            List<double[]> smoothCovariates) {
        if (rows < 1 || parametricDesign == null
                || parametricDesign.length != rows * parametricColumns
                || smoothCovariates == null
                || smoothCovariates.size() != smoothTerms.size()) {
            throw new IllegalArgumentException("prediction dimensions are invalid");
        }
        double[] beta = parametricCoefficients();
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column < parametricColumns; column++) {
                result[row] += parametricDesign[row * parametricColumns + column]
                    * beta[column];
            }
        }
        for (int term = 0; term < smoothTerms.size(); term++) {
            double[] values = smoothTerms.get(term).predict(
                smoothCovariates.get(term));
            if (values.length != rows) {
                throw new IllegalArgumentException(
                    "each smooth covariate must have one value per prediction row");
            }
            for (int row = 0; row < rows; row++) {
                result[row] += values[row];
            }
        }
        return result;
    }

    private static double[] prefix(double[] values, int length) {
        return java.util.Arrays.copyOf(values, length);
    }
}
