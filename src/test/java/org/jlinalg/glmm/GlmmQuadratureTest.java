/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

final class GlmmQuadratureTest {
    @Test
    void quadratureFitsRandomInterceptBinaryData() {
        double[] y = {0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1}; double[][] x = new double[y.length][2]; List<String> groups = List.of("a", "a", "a", "b", "b", "b", "c", "c", "c", "d", "d", "d");
        for (int row = 0; row < y.length; row++) { x[row][0] = 1.0; x[row][1] = row % 3; }
        GlmmQuadratureResult result = GlmmQuadrature.fit(y, x, groups, GlmFamilies.binomial());
        assertEquals(y.length, result.fittedMeans().length);
        assertTrue(result.converged(), result.status());
        assertTrue(result.quadratureConverged());
        assertTrue(result.gradientNorm() < 1e-5);
    }

    @Test
    void matchesIndependentRIntegralsAndLme4FitsIncludingInference() throws IOException {
        GlmmQuadratureReferenceCheck.verify(Path.of("src/test/resources/r-reference"));
    }

    @Test
    void generatedHermiteRulesIntegrateNormalMoments() {
        for (int n : new int[] {1, 2, 9, 19, 94, 257, 512}) {
            AdaptiveHermiteRule rule = AdaptiveHermiteRule.of(n);
            double mass = 0, mean = 0, second = 0, fourth = 0;
            for (int i = 0; i < n; i++) {
                double weight = Math.exp(rule.logWeights[i]) / Math.sqrt(Math.PI);
                double node = Math.sqrt(2) * rule.nodes[i];
                mass += weight; mean += weight * node; second += weight * node * node;
                fourth += weight * Math.pow(node, 4);
            }
            assertEquals(1, mass, 2e-13); assertEquals(0, mean, 2e-13);
            if (n > 1) assertEquals(1, second, 2e-12);
            if (n > 2) assertEquals(3, fourth, 2e-11);
        }
    }

    @Test
    void collapsedBinomialPreservesCombinatorialConstantAndOffsets() {
        double[][] x = {{1}};
        GlmmQuadratureEvaluation count = GlmmQuadrature.evaluate(new double[]{.5}, x, List.of("g"),
            GlmFamilies.binomial(), new double[]{100},new double[]{-.7},new double[]{.7},1,GlmmQuadratureOptions.defaults());
        assertTrue(count.converged());
        double logChoose = jdistlib.math.MathFunctions.lgammafn(101)-2*jdistlib.math.MathFunctions.lgammafn(51);
        assertEquals(-70.9414561140441435 + logChoose, count.logLikelihood(), 2e-12);
    }

    @Test
    void aSingleRuleIsLaplaceButDoesNotClaimIntegrationAccuracy() {
        double[] y = new double[100]; Arrays.fill(y, 0, 50, 1);
        double[][] x = intercepts(100);
        GlmmQuadratureEvaluation one = GlmmQuadrature.evaluate(y,x,Collections.nCopies(100,"g"),
            GlmFamilies.binomial(),new double[]{0},1,new GlmmQuadratureOptions(1,1,1e-10,100,1e-5));
        assertFalse(one.converged()); assertEquals(1,one.nodes());
        assertEquals(-100*Math.log(2)-.5*Math.log(26),one.logLikelihood(),1e-12);
        assertTrue(Double.isInfinite(one.estimatedError()));
        GlmmQuadratureEvaluation insufficient = GlmmQuadrature.evaluate(y,x,Collections.nCopies(100,"g"),
            GlmFamilies.binomial(),new double[]{0},1,new GlmmQuadratureOptions(1,5,1e-15,100,1e-5));
        assertFalse(insufficient.converged());
    }

    @Test
    void reachesExactZeroVarianceAndReportsBoundaryInference() {
        int n=40; double[] y=new double[n]; List<String> groups = new ArrayList<>();
        for (int i=0;i<n;i++) { y[i]=i%4==0?0:1; groups.add("g"+i/4); }
        GlmmQuadratureResult fit = GlmmQuadrature.fit(y,intercepts(n),groups,GlmFamilies.binomial());
        assertTrue(fit.converged(),fit.status()); assertTrue(fit.varianceBoundary());
        assertEquals(0,fit.randomVariance()); assertEquals(Math.log(3),fit.beta()[0],2e-6);
        assertEquals(Math.sqrt(1/(40*.75*.25)),fit.standardErrors()[0],2e-6);
        assertFalse(fit.jointInferenceAvailable()); assertTrue(Double.isNaN(fit.parameterCovariance()[1][1]));
        for (double mean : fit.fittedMeans()) assertEquals(.75,mean,2e-6);
        double[][] covariance=fit.fixedEffectCovariance(); covariance[0][0]=123;
        assertEquals(1/(40*.75*.25),fit.fixedEffectCovariance()[0][0],2e-6);
        double[][] scaled=intercepts(n); for(double[] row:scaled) row[0]*=1e12;
        GlmmQuadratureResult other = GlmmQuadrature.fit(y,scaled,groups,GlmFamilies.binomial());
        assertTrue(other.converged(),other.status());
        assertEquals(fit.beta()[0],other.beta()[0]*1e12,1e-8);
        assertEquals(fit.standardErrors()[0],other.standardErrors()[0]*1e12,1e-8);
    }

    @Test
    void noFalseSuccessOnSeparationOrIterationExhaustion() {
        double[] y=new double[40]; Arrays.fill(y,1); List<String> groups=new ArrayList<>();
        for(int i=0;i<40;i++) groups.add("g"+i/4);
        GlmmQuadratureResult separated=GlmmQuadrature.fit(y,intercepts(40),groups,GlmFamilies.binomial());
        assertFalse(separated.converged()); assertTrue(Double.isNaN(separated.standardErrors()[0]));
        for(int i=0;i<40;i++) y[i]=i%4==0?0:1;
        GlmmQuadratureResult limited=GlmmQuadrature.fit(y,intercepts(40),groups,GlmFamilies.binomial(),
            new GlmmQuadratureOptions(9,257,1e-10,1,1e-8));
        assertFalse(limited.converged());
    }

    @Test
    void rejectsUnsupportedFamiliesMalformedDataAndRankDeficiency() {
        double[] y={0,1,0,1}; double[][] x=intercepts(4); List<String> groups=List.of("a","a","b","b");
        for (org.jlinalg.glm.GlmFamily family : List.of(GlmFamilies.gaussian(),GlmFamilies.gamma(),
                GlmFamilies.quasiBinomial(),GlmFamilies.quasiPoisson(),GlmFamilies.negativeBinomial(2)))
            assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.fit(y,x,groups,family));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.fit(new double[]{.5,1,0,1},x,groups,GlmFamilies.binomial()));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.fit(new double[]{Double.NaN,1,0,1},x,groups,GlmFamilies.binomial()));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.fit(y,x,Collections.nCopies(4,"one"),GlmFamilies.binomial()));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.fit(y,new double[][]{{1,1},{1,1},{1,1},{1,1}},groups,GlmFamilies.binomial()));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.fit(y,x,groups,GlmFamilies.poisson(),
            new double[]{2,2,2,2},null,GlmmQuadratureOptions.defaults()));
        assertThrows(IllegalArgumentException.class,()->new GlmmQuadratureOptions(0,20,1e-8,10,1e-5));
        assertThrows(IllegalArgumentException.class,()->new GlmmQuadratureOptions(9,20,Double.NaN,10,1e-5));
    }

    private static double[][] intercepts(int n) {
        double[][] result=new double[n][1]; for(double[] row:result) row[0]=1; return result;
    }
}
