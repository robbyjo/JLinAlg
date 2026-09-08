/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.junit.jupiter.api.Test;

class LaplaceIdentificationTest {
    @Test void duplicateGroupedCovarianceRejectsBothStartsButSingleTermFits() {
        double[] y={0,0,1,0,1,1,2,1,4,5,3,6,10,12,9,11},x=new double[16];Arrays.fill(x,1);
        List<String> groups=new ArrayList<>();for(int i=0;i<16;i++)groups.add("g"+(i/4));
        var a=RandomEffectTerm.randomIntercept("a",groups);var b=RandomEffectTerm.randomIntercept("b",groups);
        var options=new GlmmLaplaceOptions(150,150,1e-8,1,1e-8,100,null);
        var fit=SparseGlmmLaplace.fit(y,x,16,1,GlmFamilies.poisson(),List.of(a),null,null,options,BackendPolicy.CPU);
        assertTrue(fit.converged());assertEquals(-30.469861767896766,fit.marginalLogLikelihood(),2e-9);
        for(double[] start:new double[][]{{1,1},{.1,2}}) {
            var controls=new GlmmLaplaceOptions(150,150,1e-8,1,1e-8,100,start);
            var ex=assertThrows(IllegalArgumentException.class,()->SparseGlmmLaplace.fit(y,x,16,1,
                GlmFamilies.poisson(),List.of(a,b),null,null,controls,BackendPolicy.CPU));
            assertTrue(ex.getMessage().contains("identifiable"));
        }
    }

    @Test void precisionRescalingAndColumnPermutationCannotConcealDuplicate() {
        double[][] a={{1,0},{1,0},{0,1},{0,1}},b={{0,2},{0,2},{2,0},{2,0}};
        var terms=List.of(RandomEffectTerm.of("a",a,List.of("x","y")),RandomEffectTerm.of("b",b,List.of("y","x")));
        var q4=new SparsePrecisionMatrix(2,new int[]{0,1,2},new int[]{0,1},new double[]{4,4});
        assertThrows(IllegalArgumentException.class,()->SparseGlmmLaplace.prepareWithPrecision(4,
            GlmFamilies.poisson(),terms,List.of(SparsePrecisionMatrix.identity(2),q4),GlmmLaplaceOptions.defaults(),BackendPolicy.CPU));
    }

    @Test void linearCombinationDependenceIsRejectedButIndependentCovariancesPass() {
        var a=RandomEffectTerm.of("a",new double[][]{{1},{1},{0},{0}},List.of("a"));
        var b=RandomEffectTerm.of("b",new double[][]{{0},{0},{1},{1}},List.of("b"));
        var sum=RandomEffectTerm.of("sum",new double[][]{{1,0},{1,0},{0,1},{0,1}},List.of("a","b"));
        try(var prepared=SparseGlmmLaplace.prepare(4,GlmFamilies.poisson(),List.of(a,b),GlmmLaplaceOptions.defaults(),BackendPolicy.CPU)) {
            assertNotNull(prepared);
        }
        assertThrows(IllegalArgumentException.class,()->SparseGlmmLaplace.prepare(4,
            GlmFamilies.poisson(),List.of(a,b,sum),GlmmLaplaceOptions.defaults(),BackendPolicy.CPU));
    }

    @Test void independentSparseCovariancesUseOnlyLinearObservationStorage() {
        int n=100_000;int[] starts=new int[n+1],first=new int[n],second=new int[n];double[] ones=new double[n];
        Arrays.fill(ones,1);for(int i=0;i<n;i++){starts[i]=i;first[i]=i/2;second[i]=i%2;}starts[n]=n;
        List<String> labels=new ArrayList<>();for(int i=0;i<n/2;i++)labels.add("g"+i);
        var a=RandomEffectTerm.ofSparseCsr("pairs",n,n/2,starts,first,ones,labels);
        var b=RandomEffectTerm.ofSparseCsr("crossed",n,2,starts,second,ones,List.of("a","b"));
        var gate=new CovarianceIdentification(n,2);
        // Identity precision solves isolate the scalable covariance actions;
        // public-API tests above cover actual sparse Cholesky integration.
        gate.add(a,v->v);gate.add(b,v->v);
        var duplicate=new CovarianceIdentification(n,2);duplicate.add(a,v->v);
        assertThrows(IllegalArgumentException.class,()->duplicate.add(a,v->v));
    }
}
