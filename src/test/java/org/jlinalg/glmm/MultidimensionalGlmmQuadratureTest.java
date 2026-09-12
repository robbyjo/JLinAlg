/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.junit.jupiter.api.Test;

final class MultidimensionalGlmmQuadratureTest {
    @Test void initialTensorGridMustFitBudgetBeforeFitOrEvaluation() {
        double[] y={0,1,0,1};double[][] x=intercepts(4);
        var effects=List.of(RandomEffectTerm.randomIntercept("g",List.of("a","a","b","b")));
        for(long budget:new long[]{1,24}) {
            var options=new MultidimensionalQuadratureOptions(5,25,budget,1e-7,200,1e-6,1e-8,1e4,20);
            var failure=assertThrows(IllegalArgumentException.class,()->
                MultidimensionalGlmmQuadrature.evaluate(y,x,effects,null,GlmFamilies.binomial(),
                    new double[]{0},new double[]{1},options,BackendPolicy.CPU));
            assertTrue(failure.getMessage().contains("maximumTotalNodes"));
            assertThrows(IllegalArgumentException.class,()->MultidimensionalGlmmQuadrature.fit(
                y,x,effects,null,GlmFamilies.binomial(),options,BackendPolicy.CPU));
        }
        var exactBudget=new MultidimensionalQuadratureOptions(5,25,25,1e-7,200,1e-6,1e-8,1e4,20);
        var result=MultidimensionalGlmmQuadrature.evaluate(y,x,effects,null,GlmFamilies.binomial(),
            new double[]{0},new double[]{1},exactBudget,BackendPolicy.CPU);
        assertEquals(25,result.nodes());assertEquals(5,result.order());
        assertTrue(Double.isFinite(result.logLikelihood()));
        assertFalse(result.converged()); // no budget remains for error estimation
    }

    @Test void overflowingTensorCountIsRejectedEvenWithMaximumLongBudget() {
        int n=30;double[] y=new double[n];double[][] x=intercepts(n);
        var groups=new java.util.ArrayList<String>();
        for(int i=0;i<n;i++){groups.add("g"+i);y[i]=i%2;}
        var effects=List.of(RandomEffectTerm.randomIntercept("g",groups));
        var options=new MultidimensionalQuadratureOptions(5,25,Long.MAX_VALUE,1e-7,200,1e-6,1e-8,1e4,20);
        assertThrows(IllegalArgumentException.class,()->MultidimensionalGlmmQuadrature.evaluate(
            y,x,effects,null,GlmFamilies.binomial(),new double[]{0},new double[]{1},options,BackendPolicy.CPU));
    }
    @Test void disjointTwoDimensionalIntegralFactorizesToScalarAghq(){
        double[] y={0,1,0,1,1,1,0,1};double[][] x=intercepts(y.length);
        double[] first=new double[y.length],second=new double[y.length];
        Arrays.fill(first,0,4,1);Arrays.fill(second,4,8,1);
        var effects=List.of(RandomEffectTerm.of("a",first,y.length,1,List.of("a")),
            RandomEffectTerm.of("b",second,y.length,1,List.of("b")));
        var options=new MultidimensionalQuadratureOptions(5,25,10000,1e-9,
            200,1e-6,1e-8,1e4,20);
        var multi=MultidimensionalGlmmQuadrature.evaluate(y,x,effects,null,
            GlmFamilies.binomial(),new double[]{0},new double[]{1,1},
            options,BackendPolicy.CPU);
        var scalar=GlmmQuadrature.evaluate(y,x,
            List.of("a","a","a","a","b","b","b","b"),
            GlmFamilies.binomial(),new double[]{0},1,
            new GlmmQuadratureOptions(5,40,1e-10,200,1e-5));
        assertTrue(multi.converged());
        assertEquals(scalar.logLikelihood(),multi.logLikelihood(),2e-8);
    }

    @Test void crossedAndPedigreePrecisionsUseTheGeneralIntegral(){
        int n=12;double[] y=new double[n];double[][] x=intercepts(n);
        String[] a=new String[n],b=new String[n];
        for(int i=0;i<n;i++){a[i]="a"+(i%2);b[i]="b"+((i/2)%2);y[i]=(i%4==0||i%5==0)?1:0;}
        var crossed=MultidimensionalGlmmQuadrature.evaluate(y,x,List.of(
            RandomEffectTerm.randomIntercept("a",Arrays.asList(a)),
            RandomEffectTerm.randomIntercept("b",Arrays.asList(b))),null,
            GlmFamilies.binomial(),new double[]{-.5},new double[]{.4,.3},
            MultidimensionalQuadratureOptions.defaults(),BackendPolicy.CPU);
        assertTrue(crossed.converged());assertTrue(Double.isFinite(crossed.logLikelihood()));

        Pedigree pedigree=Pedigree.of(List.of(PedigreeIndividual.founder("p1"),
            PedigreeIndividual.founder("p2"),new PedigreeIndividual("c","p1","p2")));
        List<String> ids=List.of("p1","p2","c","p1","p2","c","p1","p2","c","p1","p2","c");
        PedigreeRandomEffectTerm term=PedigreeRandomEffectTerm.of("pedigree",ids,pedigree);
        var related=MultidimensionalGlmmQuadrature.evaluate(y,x,
            List.of(term.randomEffect()),List.of(term.precision()),
            GlmFamilies.binomial(),new double[]{-.5},new double[]{.4},
            MultidimensionalQuadratureOptions.defaults(),BackendPolicy.CPU);
        assertTrue(related.converged());assertTrue(Double.isFinite(related.logLikelihood()));
    }

    private static double[][] intercepts(int n){double[][] x=new double[n][1];for(double[] row:x)row[0]=1;return x;}
}
