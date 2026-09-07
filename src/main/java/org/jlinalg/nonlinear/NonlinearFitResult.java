/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Fixed-effect nonlinear least-squares estimates and convergence metadata. */
public record NonlinearFitResult(
        AssociationStatistics associationStatistics,
        double[] fittedValues, double[] residuals, double objective,
        double residualVariance, int iterations, boolean converged,
        BackendProvenance backend) {
    public NonlinearFitResult {
        fittedValues = fittedValues.clone(); residuals = residuals.clone();
    }
    public double[] beta() { return associationStatistics.beta(); }
    public double[] standardErrors() { return associationStatistics.standardErrors(); }
    @Override public double[] fittedValues() { return fittedValues.clone(); }
    @Override public double[] residuals() { return residuals.clone(); }
}
