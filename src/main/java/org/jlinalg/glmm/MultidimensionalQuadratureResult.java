/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Maximum marginal-likelihood fit from multidimensional adaptive quadrature. */
public record MultidimensionalQuadratureResult(String family,double[] fixedEffects,
        double[] varianceComponents,double logLikelihood,int quadratureOrder,
        long quadratureNodes,double quadratureError,boolean quadratureConverged,
        int evaluations,boolean converged,String status) {
    public MultidimensionalQuadratureResult{
        fixedEffects=fixedEffects.clone();varianceComponents=varianceComponents.clone();
    }
    @Override public double[] fixedEffects(){return fixedEffects.clone();}
    @Override public double[] varianceComponents(){return varianceComponents.clone();}
}
