/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.inference;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

class CoreInferenceAuditTest {
    @Test void nonfiniteCoefficientsCannotProduceSignificance() {
        for(double beta:new double[]{Double.NaN,Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY})
            for(double se:new double[]{0,1}){
                var fit=AssociationStatistics.normal(new double[]{beta},new double[]{se});
                assertTrue(Double.isNaN(fit.pValues()[0]));
                assertTrue(Double.isNaN(fit.negativeLog10PValues()[0]));
            }
        assertEquals(0,AssociationStatistics.normal(new double[]{1},new double[]{0}).pValues()[0]);
    }
    @Test void contrastRejectsNonfiniteOrAsymmetricCovariance() {
        assertThrows(IllegalArgumentException.class,()->LinearHypothesis.test(new double[]{1,2},
            new double[]{1,0.5,0,1},new double[][]{{1,0},{0,1}},20,BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class,()->LinearHypothesis.test(new double[]{Double.NaN},
            new double[]{1},new double[][]{{1}},20,BackendPolicy.CPU));
        var fit=LinearHypothesis.test(new double[]{1,2},new double[]{2,0.5,0.5,1},
            new double[][]{{1,0},{0,1}},20,BackendPolicy.CPU);
        assertEquals(2,fit.statistic(),1e-14);
    }
    @Test void modelDimensionsDoNotWrapIntegerProducts() {
        assertThrows(IllegalArgumentException.class,()->MatrixOps.validateModelData(
            new double[65536],new double[0],65536,65536));
    }
    @Test void bootstrapSummariesAreScaleSafeAndValidateAccounting() {
        var result=new GaussianBootstrapResult(3,1,.95,new double[]{2e200},java.util.List.of(),
            new double[0],new double[][]{{1e200},{2e200},{3e200}},new double[3][0],java.util.List.of());
        var summary=result.fixedEffectSummaries().get(0);
        assertEquals(2,summary.bootstrapMean()/1e200,1e-14);
        assertEquals(1,summary.standardError()/1e200,1e-14);
        assertThrows(IllegalArgumentException.class,()->new GaussianBootstrapResult(4,1,.95,
            new double[]{1},java.util.List.of(),new double[0],new double[][]{{1},{2},{3}},
            new double[3][0],java.util.List.of()));
        assertThrows(IllegalArgumentException.class,()->new GaussianBootstrapResult(2,1,.95,
            new double[]{1},java.util.List.of(),new double[0],new double[][]{{1},{Double.NaN}},
            new double[2][0],java.util.List.of()));
    }
    @Test void constantMaximumBootstrapDrawsRemainExactlyConstant() {
        double[][] draws=new double[9][1];for(double[] row:draws)row[0]=Double.MAX_VALUE;
        var fit=new GaussianBootstrapResult(9,1,.95,new double[]{Double.MAX_VALUE},java.util.List.of(),
            new double[0],draws,new double[9][0],java.util.List.of());
        var summary=fit.fixedEffectSummaries().get(0);
        assertEquals(Double.MAX_VALUE,summary.bootstrapMean());
        assertEquals(0,summary.bias());assertEquals(0,summary.standardError());
    }
    @Test void waldContrastsAreUnitInvariantAndRequireFullPsdCovariance() {
        for(double scale:new double[]{1e-200,1,1e200}){
            var fit=LinearHypothesis.test(new double[]{1,2},new double[]{1,0,0,1},
                new double[][]{{scale,0},{0,scale}},20,BackendPolicy.CPU);
            assertEquals(2.5,fit.statistic(),1e-14);
        }
        assertThrows(IllegalArgumentException.class,()->LinearHypothesis.test(new double[]{1,2},
            new double[]{1,2,2,1},new double[][]{{1,0}},20,BackendPolicy.CPU));
        var singular=LinearHypothesis.test(new double[]{1,2},new double[]{1,1,1,1},
            new double[][]{{1,1}},20,BackendPolicy.CPU);
        assertEquals(2.25,singular.statistic(),1e-14);
    }
}
