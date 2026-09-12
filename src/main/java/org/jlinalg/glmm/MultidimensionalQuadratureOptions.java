/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Controls for low-dimensional tensor adaptive Gauss-Hermite integration. */
public record MultidimensionalQuadratureOptions(int initialOrder,int maximumOrder,
        long maximumTotalNodes,double quadratureTolerance,
        int maximumOuterEvaluations,double relativeTolerance,
        double minimumVariance,double maximumVariance,
        double maximumAbsoluteCoefficient, IntegrationMethod integrationMethod) {
    /** Tensor rule or reproducible randomly shifted Halton importance integration. */
    public enum IntegrationMethod { ADAPTIVE_HERMITE, SHIFTED_HALTON }
    public MultidimensionalQuadratureOptions(int initialOrder,int maximumOrder,
            long maximumTotalNodes,double quadratureTolerance,int maximumOuterEvaluations,
            double relativeTolerance,double minimumVariance,double maximumVariance,
            double maximumAbsoluteCoefficient) {
        this(initialOrder,maximumOrder,maximumTotalNodes,quadratureTolerance,
            maximumOuterEvaluations,relativeTolerance,minimumVariance,maximumVariance,
            maximumAbsoluteCoefficient,IntegrationMethod.ADAPTIVE_HERMITE);
    }
    public MultidimensionalQuadratureOptions{
        if(initialOrder<1||maximumOrder<initialOrder||maximumTotalNodes<1
                ||!(quadratureTolerance>0)||maximumOuterEvaluations<20
                ||!(relativeTolerance>0)||!(minimumVariance>0)
                ||!(maximumVariance>minimumVariance)
                ||!(maximumAbsoluteCoefficient>0)||integrationMethod==null
                ||!Double.isFinite(quadratureTolerance)||!Double.isFinite(relativeTolerance)
                ||!Double.isFinite(maximumVariance)||!Double.isFinite(maximumAbsoluteCoefficient))
            throw new IllegalArgumentException("invalid multidimensional quadrature controls");
    }
    public static MultidimensionalQuadratureOptions defaults(){
        return new MultidimensionalQuadratureOptions(5,25,2_000_000,
            1e-7,500,1e-6,1e-8,1e8,30);
    }
    /** Orders are points per replicate; eight fixed independent random shifts.
     * Tolerance uses a three-standard-error estimate, not a certified bound. */
    public static MultidimensionalQuadratureOptions nonTensorDefaults(){
        return new MultidimensionalQuadratureOptions(1024,65536,524288,
            1e-4,500,1e-5,1e-8,1e8,30,IntegrationMethod.SHIFTED_HALTON);
    }
}
