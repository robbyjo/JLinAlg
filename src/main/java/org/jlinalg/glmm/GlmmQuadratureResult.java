/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Independent random-intercept adaptive-quadrature fit.
 * Fitted means average over a NEW group's random intercept (not posterior modes).
 * Binomial means are probabilities, including when trials are supplied; a mean
 * whose own integration budget is exhausted is NaN.
 * Joint covariance is ordered fixed coefficients, random VARIANCE, and includes
 * nuisance-variance uncertainty at an interior optimum. At zero variance its
 * variance row/column is NaN; fixed covariance conditions on variance zero.
 * Nonconverged fits have NaN covariance. Wald inference is asymptotic.
 * Gradient norm is the projected score in internally scaled fixed-coefficient
 * and log1p(variance / within-group variance scale) coordinates. Reported
 * parameter covariance is transformed back to coefficients and variance. */
public record GlmmQuadratureResult(String family, double[] fixedEffects,
                                   double randomStandardDeviation, double[] fittedMeans,
                                   double logLikelihood, int nodes, int evaluations,
                                   boolean converged, double[][] fixedEffectCovariance,
                                   double[][] parameterCovariance, double gradientNorm,
                                   double quadratureError, boolean quadratureConverged,
                                   boolean jointInferenceAvailable, boolean varianceBoundary,
                                   String status) {
    public GlmmQuadratureResult {
        fixedEffects = fixedEffects.clone(); fittedMeans = fittedMeans.clone();
        fixedEffectCovariance = copy(fixedEffectCovariance); parameterCovariance = copy(parameterCovariance);
    }
    /** Source-compatible constructor for legacy externally constructed results;
     * no inference or integration accuracy is implied by absent diagnostics. */
    public GlmmQuadratureResult(String family, double[] fixedEffects, double randomStandardDeviation,
            double[] fittedMeans, double logLikelihood, int nodes, int evaluations, boolean converged) {
        this(family, fixedEffects, randomStandardDeviation, fittedMeans, logLikelihood, nodes, evaluations,
            converged, unavailable(fixedEffects.length), unavailable(fixedEffects.length + 1),
            Double.NaN, Double.NaN, false, false, randomStandardDeviation == 0, "legacy result without diagnostics");
    }
    public double[] fixedEffects() { return fixedEffects.clone(); }
    public double[] beta() { return fixedEffects(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
    public double[][] fixedEffectCovariance() { return copy(fixedEffectCovariance); }
    public double[][] parameterCovariance() { return copy(parameterCovariance); }
    public double randomVariance() { return randomStandardDeviation * randomStandardDeviation; }
    public double[] standardErrors() {
        double[] result = new double[fixedEffects.length];
        for (int i = 0; i < result.length; i++) result[i] = Math.sqrt(fixedEffectCovariance[i][i]);
        return result;
    }
    public double[] waldZ() {
        double[] result = standardErrors();
        for (int i = 0; i < result.length; i++) result[i] = fixedEffects[i] / result[i];
        return result;
    }
    public double[] pValues() {
        double[] result = waldZ();
        for (int i = 0; i < result.length; i++)
            result[i] = 2 * jdistlib.Normal.cumulative(-Math.abs(result[i]), 0, 1, true, false);
        return result;
    }
    private static double[][] copy(double[][] a) {
        double[][] result = new double[a.length][];
        for (int i = 0; i < a.length; i++) result[i] = a[i].clone();
        return result;
    }
    private static double[][] unavailable(int n) {
        double[][] result = new double[n][n];
        for (double[] row : result) java.util.Arrays.fill(row, Double.NaN);
        return result;
    }
}
