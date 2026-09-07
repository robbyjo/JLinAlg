/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.mixed.SparseLinearMixedModelResult;

/** Nonlinear fixed effects with sparse ordinary or pedigree random effects. */
public record NonlinearMixedFitResult(
        double[] parameters, double[] fittedValues, double[] residuals,
        double objective, int iterations, boolean converged,
        SparseLinearMixedModelResult linearizedModel,
        BackendProvenance backend) {
    public NonlinearMixedFitResult {
        parameters = parameters.clone(); fittedValues = fittedValues.clone(); residuals = residuals.clone();
    }
    @Override public double[] parameters() { return parameters.clone(); }
    @Override public double[] fittedValues() { return fittedValues.clone(); }
    @Override public double[] residuals() { return residuals.clone(); }
}
