/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** One low-dimensional adaptive Gaussian random-effect integral. */
public record MultidimensionalQuadratureEvaluation(double logLikelihood,
        int order,long nodes,double estimatedError,boolean converged,
        double[] posteriorMode) {
    public MultidimensionalQuadratureEvaluation{posteriorMode=posteriorMode.clone();}
    @Override public double[] posteriorMode(){return posteriorMode.clone();}
}
