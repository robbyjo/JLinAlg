/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** FIML missing-pattern estimate together with the fitted SEM. */
public record SemFimlResult(SemFitResult fit, double[] means, double[] covariance,
                            int patternCount, int iterations, boolean converged) {
    public SemFimlResult { means = means.clone(); covariance = covariance.clone(); }
    public double[] means() { return means.clone(); }
    public double[] covariance() { return covariance.clone(); }
}
