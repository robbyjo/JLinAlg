/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** Full-covariance GLS meta-regression result. */
public record MetaMultilevelRegressionResult(java.util.List<String> coefficientNames,
        double[] coefficients, double[] covariance, double q, double degreesOfFreedom,
        double logLikelihood, boolean converged) {
    public MetaMultilevelRegressionResult { coefficientNames = java.util.List.copyOf(coefficientNames); coefficients = coefficients.clone(); covariance = covariance.clone(); }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] beta() { return coefficients(); }
    public double[] covariance() { return covariance.clone(); }
    public double[] standardErrors() { double[] result = new double[coefficients.length]; for (int i = 0; i < result.length; i++) result[i] = Math.sqrt(Math.max(0.0, covariance[i * result.length + i])); return result; }
}
