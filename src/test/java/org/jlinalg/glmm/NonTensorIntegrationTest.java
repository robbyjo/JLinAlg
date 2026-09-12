/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;
import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.mixed.RandomEffectTerm;
import org.junit.jupiter.api.Test;
final class NonTensorIntegrationTest {
    @Test void eightDimensionalIntegralAvoidsTensorGrowth() {
        int n=8;double[] y=new double[n];double[][] x=new double[n][1];java.util.ArrayList<String> groups=new java.util.ArrayList<>();
        for(int i=0;i<n;i++){y[i]=i%2;x[i][0]=1;groups.add("g"+i);}
        var options=new MultidimensionalQuadratureOptions(1024,65536,524288,2e-3,100,1e-5,1e-8,1e4,30,MultidimensionalQuadratureOptions.IntegrationMethod.SHIFTED_HALTON);
        var result=MultidimensionalGlmmQuadrature.evaluate(y,x,List.of(RandomEffectTerm.randomIntercept("g",groups)),null,GlmFamilies.binomial(),new double[]{0},new double[]{1},options,BackendPolicy.CPU);
        assertTrue(result.converged(),result.toString());assertEquals(-8*Math.log(2),result.logLikelihood(),2e-3);
        assertTrue(result.nodes()<=524288);
    }
    @Test void importanceNormalizerIntegratesGaussianExactly() {
        var estimate=ShiftedHalton.integrate(3,100,new double[3],new double[]{1,0,0,0,1,0,0,0,1},
            x->-.5*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]),0);
        assertEquals(0,estimate.logIntegral(),1e-14);assertEquals(0,estimate.error(),1e-14);
    }
    @Test void nonTensorLikelihoodMatchesIndependentScalarIntegralAndRepeatsExactly() {
        double[] y={0,1,0,1,1,1,0,1};double[][] x=new double[y.length][1];for(double[] row:x)row[0]=1;
        var groups=List.of("a","a","a","a","b","b","b","b");
        var effect=List.of(RandomEffectTerm.randomIntercept("g",groups));
        var options=MultidimensionalQuadratureOptions.nonTensorDefaults();
        var value=MultidimensionalGlmmQuadrature.evaluate(y,x,effect,null,GlmFamilies.binomial(),new double[]{0},new double[]{1},options,BackendPolicy.CPU);
        var scalar=GlmmQuadrature.evaluate(y,x,groups,GlmFamilies.binomial(),new double[]{0},1,new GlmmQuadratureOptions(5,40,1e-10,200,1e-5));
        assertTrue(value.converged(),value.toString());assertEquals(scalar.logLikelihood(),value.logLikelihood(),2e-4);
        var again=MultidimensionalGlmmQuadrature.evaluate(y,x,effect,null,GlmFamilies.binomial(),new double[]{0},new double[]{1},options,BackendPolicy.CPU);
        assertEquals(value.logLikelihood(),again.logLikelihood());assertEquals(value.estimatedError(),again.estimatedError());
        var tiny=new MultidimensionalQuadratureOptions(16,16,128,1e-12,100,1e-6,1e-8,1e4,30,MultidimensionalQuadratureOptions.IntegrationMethod.SHIFTED_HALTON);
        assertFalse(MultidimensionalGlmmQuadrature.evaluate(y,x,effect,null,GlmFamilies.binomial(),new double[]{0},new double[]{1},tiny,BackendPolicy.CPU).converged());
    }
}
