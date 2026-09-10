/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Controls for low-dimensional tensor adaptive Gauss-Hermite integration. */
public record MultidimensionalQuadratureOptions(int initialOrder,int maximumOrder,
        long maximumTotalNodes,double quadratureTolerance,
        int maximumOuterEvaluations,double relativeTolerance,
        double minimumVariance,double maximumVariance,
        double maximumAbsoluteCoefficient) {
    public MultidimensionalQuadratureOptions{
        if(initialOrder<1||maximumOrder<initialOrder||maximumTotalNodes<1
                ||!(quadratureTolerance>0)||maximumOuterEvaluations<20
                ||!(relativeTolerance>0)||!(minimumVariance>0)
                ||!(maximumVariance>minimumVariance)
                ||!(maximumAbsoluteCoefficient>0))
            throw new IllegalArgumentException("invalid multidimensional quadrature controls");
    }
    public static MultidimensionalQuadratureOptions defaults(){
        return new MultidimensionalQuadratureOptions(5,25,2_000_000,
            1e-7,500,1e-6,1e-8,1e8,30);
    }
}
