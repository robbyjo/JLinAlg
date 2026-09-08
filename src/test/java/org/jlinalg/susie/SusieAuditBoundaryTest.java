/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class SusieAuditBoundaryTest {
    private final SusieOptions fixed=new SusieOptions(1,20,1e-8,.2,false,.95,0);
    @Test void tinyPositivePriorRetainsRepresentableConjugateMean() {
        for(double sign:new double[]{-1,1}) {
            var fit=Susie.fitSufficientStatistics(new double[]{1e-200},new double[]{sign*1e-100},1,2,null,
                new SusieOptions(1,20,1e-8,1e-200,false,.95,0),BackendPolicy.CPU);
            assertTrue(fit.converged());
            // Posterior variance V/(1+V*XtX)=1e-200; posterior mean=variance*Xty.
            assertEquals(sign,fit.posteriorMean()[0]/1e-300,1e-13);
            assertEquals(1,fit.pip()[0]);
        }
    }
    @Test void extremeZHasAnalyticSaturatedEffectAndNonzeroTinyPip() {
        var fit=Susie.fitSummary(new double[]{1e200,0},new double[][]{{1,0},{0,1}},100,null,fixed,BackendPolicy.CPU);
        assertTrue(fit.converged());
        double shrink=.2/(.2+1.0/99);
        double weak=1/(1+Math.exp(.5*99*shrink));
        assertEquals(shrink,fit.posteriorMean()[0],1e-15);
        assertEquals(weak,fit.pip()[1],weak*2e-14);
        assertTrue(fit.pip()[1]>0);
    }
    @Test void iterationLimitReturnsTheVarianceThatGeneratedItsPosterior() {
        var fit=Susie.fitSufficientStatistics(new double[]{10},new double[]{5},10,11,null,
            new SusieOptions(1,1,1e-8,.2,true,.95,0),BackendPolicy.CPU);
        assertFalse(fit.converged()); assertEquals(1,fit.iterations());
        assertEquals(1,fit.residualVariance());
        assertEquals(.2/(.2+fit.residualVariance()/10)*.5,fit.posteriorMean()[0],1e-15);
    }
    @Test void rejectsIndefiniteLdButAcceptsSingularLd() {
        assertThrows(IllegalArgumentException.class,()->Susie.fitSummary(new double[]{2,3,4},
            new double[][]{{1,.9,.9},{.9,1,-.9},{.9,-.9,1}},100,null,fixed,BackendPolicy.CPU));
        var fit=Susie.fitSummary(new double[]{3,3},new double[][]{{1,1},{1,1}},100,null,fixed,BackendPolicy.CPU);
        assertTrue(fit.converged()); assertArrayEquals(new double[]{.5,.5},fit.pip(),1e-15);
    }
    @Test void rejectsInvalidScalesDimensionsAndDuplicateNames() {
        for(double n:new double[]{2,10.5,Double.POSITIVE_INFINITY,3e10})
            assertThrows(IllegalArgumentException.class,()->Susie.fitSummary(new double[]{1},new double[][]{{1}},n,null,fixed,BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class,()->Susie.fitSummary(new double[]{1,2},new double[][]{{1},{1}},100,null,fixed,BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class,()->Susie.fitSummary(new double[]{1,2},new double[][]{{1,0},{0,1}},100,List.of("a","a"),fixed,BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class,()->new SusieOptions(1,2,1e-6,Double.POSITIVE_INFINITY,false,.95,0));
        assertThrows(IllegalArgumentException.class,()->Susie.fitSufficientStatistics(new double[]{1},new double[]{Double.NaN},10,11,null,fixed,BackendPolicy.CPU));
        assertThrows(IllegalArgumentException.class,()->Susie.fitSufficientStatistics(new double[]{1},new double[]{100},1,11,null,fixed,BackendPolicy.CPU));
    }
}
