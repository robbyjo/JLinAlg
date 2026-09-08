/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.reml;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.junit.jupiter.api.Test;
class RemainingRemlAccuracyTest {
    @Test void unbalancedMlAndRemlMatchIndependentR() throws Exception {
        Path root=Path.of("src/test/resources/r-reference/remaining-mixed");var p=new Properties();
        try(var reader=Files.newBufferedReader(root.resolve("reference.properties"))){p.load(reader);}
        var lines=Files.readAllLines(root.resolve("reml.tsv"));int n=lines.size()-1;double[] y=new double[n];double[][] x=new double[n][2];List<String> groups=new ArrayList<>();
        for(int i=0;i<n;i++){var a=lines.get(i+1).split("\t");y[i]=Double.parseDouble(a[0]);x[i][0]=1;x[i][1]=Double.parseDouble(a[1]);groups.add(a[2]);}
        for(var method:VarianceEstimation.values()) {
            var fit=Reml.fit(y,x,List.of(VarianceComponent.randomIntercept("g",groups),VarianceComponent.identity("residual",n)),
                RemlOptions.builder().varianceEstimation(method).build(),BackendPolicy.CPU);
            assertTrue(fit.converged(),fit.convergenceMessage());String prefix="reml."+method;
            assertArrayEquals(numbers(p,prefix+".variance"),fit.varianceComponents(),2e-6);
            assertArrayEquals(numbers(p,prefix+".beta"),fit.beta(),2e-6);
            assertArrayEquals(numbers(p,prefix+".covariance"),fit.fixedEffectCovariance(),2e-6);
            assertEquals(numbers(p,prefix+".LL")[0],fit.restrictedLogLikelihood(),2e-8);
        }
    }
    @Test void redundantBasesCannotAcquireIdentificationThroughRidge() {
        double[] y={-1,1,-2,2};double[][] x={{1},{1},{1},{1}};
        for(var method:List.of(DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION,DegreesOfFreedomMethod.SATTERTHWAITE,DegreesOfFreedomMethod.KENWARD_ROGER))
            assertThrows(IllegalArgumentException.class,()->Reml.fit(y,x,List.of(VarianceComponent.identity("a",4),VarianceComponent.identity("b",4)),
                RemlOptions.builder().degreesOfFreedomMethod(method).build(),BackendPolicy.CPU));
    }
    @Test void knownPositiveCovarianceCannotMaskAnIndefiniteVarianceBasis() {
        var bad=new VarianceComponent("invalid",4,new double[]{-1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1});
        double[] known={10,0,0,0,0,10,0,0,0,0,10,0,0,0,0,10};
        assertThrows(IllegalArgumentException.class,()->Reml.fitWithKnownCovariance(new double[]{-1,1,-2,2},new double[]{1,1,1,1},4,1,List.of(bad),known,
            RemlOptions.builder().initialVariances(.1).build(),BackendPolicy.CPU));
    }
    @Test void remlCannotIdentifyCovarianceInsideTheFixedMeanSpace() {
        assertThrows(IllegalArgumentException.class,()->Reml.fit(new double[]{-1,1,-2,2},new double[][]{{1},{1},{1},{1}},
            List.of(VarianceComponent.randomIntercept("oneGroup",List.of(1,1,1,1)),VarianceComponent.identity("residual",4)),
            RemlOptions.defaults(),BackendPolicy.CPU));
    }
    private static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
