/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import jdistlib.Normal;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.internal.MatrixOps;

/** Estimated marginal means and treatment/visit contrasts for GEE reports. */
public final class ClinicalGee {
    private ClinicalGee() { }

    /** Computes one adjusted mean for each supplied reference-grid design row. */
    public static MarginalMeanEstimate[] marginalMeans(
            GeeResult fit,
            double[][] referenceGrid,
            GlmFamily family) {
        return marginalMeans(fit, referenceGrid, family, 0.95);
    }

    /** Computes adjusted means and delta-method uncertainty. */
    public static MarginalMeanEstimate[] marginalMeans(
            GeeResult fit,
            double[][] referenceGrid,
            GlmFamily family,
            double confidenceLevel) {
        if (fit == null || family == null || referenceGrid == null
                || referenceGrid.length == 0
                || !(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException("marginal-means inputs are invalid");
        }
        double[] beta = fit.coefficients();
        double[] covariance = fit.covariance();
        double critical = Normal.quantile(0.5 + confidenceLevel / 2.0,
            0.0, 1.0, true, false);
        MarginalMeanEstimate[] result = new MarginalMeanEstimate[referenceGrid.length];
        for (int row = 0; row < referenceGrid.length; row++) {
            double[] design = validateRow(referenceGrid[row], beta.length);
            double eta = dot(design, beta);
            double variance = quadratic(design, covariance);
            double linkStandardError = Math.sqrt(Math.max(0.0, variance));
            double mean = family.inverseLink(eta);
            double responseStandardError = Math.abs(family.meanDerivative(eta))
                * linkStandardError;
            double lower = family.inverseLink(eta - critical * linkStandardError);
            double upper = family.inverseLink(eta + critical * linkStandardError);
            result[row] = new MarginalMeanEstimate(eta, linkStandardError,
                mean, responseStandardError, Math.min(lower, upper),
                Math.max(lower, upper));
        }
        return result;
    }

    /** Computes a link-scale contrast between two reference-grid rows. */
    public static MarginalContrast contrast(
            GeeResult fit,
            double[] first,
            double[] second) {
        return contrast(fit, first, second, 0.95);
    }

    /** Computes a link-scale contrast, Wald test, and exponentiated effect. */
    public static MarginalContrast contrast(
            GeeResult fit,
            double[] first,
            double[] second,
            double confidenceLevel) {
        if (fit == null || !(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException("contrast inputs are invalid");
        }
        double[] beta = fit.coefficients();
        double[] left = validateRow(first, beta.length);
        double[] right = validateRow(second, beta.length);
        double[] contrast = new double[beta.length];
        for (int index = 0; index < contrast.length; index++) {
            contrast[index] = left[index] - right[index];
        }
        double estimate = dot(contrast, beta);
        double standardError = Math.sqrt(Math.max(0.0,
            quadratic(contrast, fit.covariance())));
        double statistic = estimate / standardError;
        double pValue = Double.isFinite(statistic)
            ? 2.0 * Normal.cumulative(Math.abs(statistic),
                0.0, 1.0, false, false)
            : Double.NaN;
        double critical = Normal.quantile(0.5 + confidenceLevel / 2.0,
            0.0, 1.0, true, false);
        double lower = estimate - critical * standardError;
        double upper = estimate + critical * standardError;
        return new MarginalContrast(estimate, standardError, statistic, pValue,
            lower, upper, safeExp(estimate), safeExp(lower), safeExp(upper));
    }

    private static double[] validateRow(double[] row, int columns) {
        if (row == null || row.length != columns) {
            throw new IllegalArgumentException(
                "reference-grid columns must equal coefficient count");
        }
        return MatrixOps.finiteCopy(row, "reference-grid row");
    }

    private static double dot(double[] left, double[] right) {
        double result = 0.0;
        for (int index = 0; index < left.length; index++) {
            result += left[index] * right[index];
        }
        return result;
    }

    private static double quadratic(double[] vector, double[] matrix) {
        int dimension = vector.length;
        double result = 0.0;
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                result += vector[row] * matrix[row * dimension + column]
                    * vector[column];
            }
        }
        return result;
    }

    private static double safeExp(double value) {
        return Math.exp(Math.max(-700.0, Math.min(700.0, value)));
    }
}
