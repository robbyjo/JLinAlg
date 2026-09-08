/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Properties;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

/** Independent review regressions; frozen expectations come from base R. */
final class GlmmQuadratureStressTest {
    private static final Path ROOT=Path.of("src/test/resources/r-reference");
    private static final GlmmQuadratureOptions OPTIONS=GlmmQuadratureOptions.defaults();
    private static Properties reference() throws IOException {
        Properties p=new Properties();
        try(Reader reader=Files.newBufferedReader(ROOT.resolve("quadrature-stress.properties"))) {p.load(reader);}
        return p;
    }
    private static double number(Properties p,String key) {return Double.parseDouble(p.getProperty(key));}

    @Test void concentratedOptimumMustNotBeReportedAsZeroVariance() throws IOException {
        Properties p=reference(); double n=number(p,"boundary.n"), k=number(p,"boundary.k");
        double[] y={(n/2-k)/n,(n/2+k)/n},trials={n,n}; double[][] x={{1},{1}};
        GlmmQuadratureResult fit=GlmmQuadrature.fit(y,x,List.of("a","b"),GlmFamilies.binomial(),trials,null,OPTIONS);
        assertTrue(fit.converged(),fit.status()); assertFalse(fit.varianceBoundary()); assertTrue(fit.jointInferenceAvailable());
        assertEquals(0,fit.beta()[0],1e-8);
        assertEquals(number(p,"boundary.variance"),fit.randomVariance(),3e-11);
        assertEquals(number(p,"boundary.logLikelihood"),fit.logLikelihood(),2e-9);
        assertTrue(fit.logLikelihood()>number(p,"boundary.zeroLogLikelihood")+.3);
        assertEquals(number(p,"boundary.betaSE"),fit.standardErrors()[0],2e-8);
        assertEquals(number(p,"boundary.varianceCovariance"),fit.parameterCovariance()[1][1],2e-15);
        n=number(p,"shallow.n");k=number(p,"shallow.k");
        fit=GlmmQuadrature.fit(new double[]{(n/2-k)/n,(n/2+k)/n},x,List.of("a","b"),
            GlmFamilies.binomial(),new double[]{n,n},null,OPTIONS);
        assertTrue(fit.converged(),fit.status());assertFalse(fit.varianceBoundary());
        assertEquals(number(p,"shallow.variance"),fit.randomVariance(),3e-11);
        assertEquals(number(p,"shallow.logLikelihood"),fit.logLikelihood(),2e-9);
        assertTrue(fit.logLikelihood()>number(p,"shallow.zeroLogLikelihood")+.004);
    }

    @Test void concentratedVarianceScaleGridRetainsInteriorSolutions() {
        for(int n:new int[]{10000,100000,300000,1000000,3000000}) {
            double k=Math.rint(Math.sqrt(n/2.0)); double[] y={(n/2.-k)/n,(n/2.+k)/n},trials={n,n};
            double[][] x={{1},{1}}; var groups=List.of("a","b");
            var fit=GlmmQuadrature.fit(y,x,groups,GlmFamilies.binomial(),trials,null,OPTIONS);
            var zero=GlmmQuadrature.evaluate(y,x,groups,GlmFamilies.binomial(),trials,null,new double[]{0},0,OPTIONS);
            assertTrue(fit.converged(),"n="+n+" "+fit.status());
            assertFalse(fit.varianceBoundary(),"n="+n);
            assertTrue(fit.logLikelihood()>zero.logLikelihood()+.3,"n="+n);
        }
    }

    @Test void largeCountProbabilitiesRemainAccurateAndNonpositive() throws IOException {
        for(String line:Files.readAllLines(ROOT.resolve("quadrature-large-counts.tsv")).subList(1,10)) {
            String[] a=line.split("\t"); double n=Double.parseDouble(a[0]),s=Math.floor(n/2);
            var bin=GlmmQuadrature.evaluate(new double[]{s/n},new double[][]{{1}},List.of("g"),GlmFamilies.binomial(),
                new double[]{n},null,new double[]{0},0,OPTIONS);
            var poisson=GlmmQuadrature.evaluate(new double[]{n},new double[][]{{1}},List.of("g"),GlmFamilies.poisson(),
                new double[]{Math.log(n)},0,OPTIONS);
            assertTrue(bin.converged());assertTrue(poisson.converged());
            assertTrue(bin.logLikelihood()<=0);assertTrue(poisson.logLikelihood()<=0);
            assertEquals(Double.parseDouble(a[1]),bin.logLikelihood(),2e-11,"binomial n="+n);
            assertEquals(Double.parseDouble(a[2]),poisson.logLikelihood(),2e-11,"Poisson n="+n);
        }
    }

    @Test void validLargeDenominatorProportionsSurviveRoundTripButFractionsAreRejected() throws IOException {
        Properties p=reference();double n=number(p,"proportion.trials"),s=number(p,"proportion.successes");
        double y=s/n;
        assertTrue(Math.abs(y*n-s)>1e-8,"counterexample must exceed old absolute tolerance");
        var result=GlmmQuadrature.evaluate(new double[]{y},new double[][]{{1}},List.of("g"),GlmFamilies.binomial(),
            new double[]{n},null,new double[]{Math.log(s/(n-s))},0,OPTIONS);
        assertTrue(result.converged());assertEquals(number(p,"proportion.logLikelihood"),result.logLikelihood(),2e-11);
        for(double invalid:new double[]{(s+.1)/n,(s+.25)/n})
            assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.evaluate(new double[]{invalid},new double[][]{{1}},List.of("g"),
                GlmFamilies.binomial(),new double[]{n},null,new double[]{0},0,OPTIONS));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.evaluate(new double[]{1e-9},new double[][]{{1}},List.of("g"),
            GlmFamilies.binomial(),new double[]{0},0,OPTIONS));
        assertThrows(IllegalArgumentException.class,()->GlmmQuadrature.evaluate(new double[]{1+1e-9},new double[][]{{1}},List.of("g"),
            GlmFamilies.poisson(),new double[]{0},0,OPTIONS));
    }

    @Test void bracketedModesConvergeAcrossLargeOffsets() throws IOException {
        Properties p=reference();
        for(int offset:new int[]{100,300,500,1000}) {
            var result=GlmmQuadrature.evaluate(new double[]{1},new double[][]{{1}},List.of("g"),GlmFamilies.poisson(),
                null,new double[]{offset},new double[]{0},1,OPTIONS);
            assertTrue(result.converged(),"offset="+offset+" "+result);
            assertEquals(number(p,"offset."+offset),result.logLikelihood(),2e-9,"offset="+offset);
        }
    }
}
