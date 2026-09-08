/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** One-dimensional random-intercept Gaussian-quadrature GLMM result. */
public record GlmmQuadratureResult(String family, double[] fixedEffects,
                                   double randomStandardDeviation, double[] fittedMeans,
                                   double logLikelihood, int nodes, int evaluations,
                                   boolean converged) {
    public GlmmQuadratureResult { fixedEffects = fixedEffects.clone(); fittedMeans = fittedMeans.clone(); }
    public double[] fixedEffects() { return fixedEffects.clone(); }
    public double[] beta() { return fixedEffects(); }
    public double[] fittedMeans() { return fittedMeans.clone(); }
}
