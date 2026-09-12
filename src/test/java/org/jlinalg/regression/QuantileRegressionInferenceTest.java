/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;
import org.junit.jupiter.api.Test;

class QuantileRegressionInferenceTest {
    @Test void densityAndKernelCovarianceMatchIndependentRqAndMatrixReference() throws Exception {
        var d=ExactQuantileTest.data("regular");Properties ref=new Properties();
        try(var in=Files.newInputStream(Path.of("src/test/resources/feasible-extensions/quantile.properties"))){ref.load(in);}
        double tau=numbers(ref,"tau")[0];
        var fit=QuantileRegressionInference.fitExact(d.y(),d.x(),tau,numbers(ref,"densities"));
        assertArrayEquals(numbers(ref,"estimates"),fit.fit().coefficients(),2e-6);
        assertArrayEquals(numbers(ref,"covariance"),fit.coefficientCovariance(),1e-12);
        var kernel=QuantileRegressionInference.fitExactIidKernel(d.y(),d.x(),tau,numbers(ref,"bandwidth")[0]);
        assertEquals(numbers(ref,"kernelDensity")[0],kernel.conditionalDensities()[0],2e-8);
        assertArrayEquals(numbers(ref,"kernelCovariance"),kernel.coefficientCovariance(),2e-8);
    }
    @Test void normalMedianHasAnalyticVarianceAndCovarianceRespectsUnits() throws Exception {
        int n=101;double[] y=new double[n],density=new double[n];double[][] x=new double[n][1];
        Arrays.fill(density,1/Math.sqrt(2*Math.PI));
        for(int i=0;i<n;i++){y[i]=jdistlib.Normal.quantile((i+.5)/n,0,1,true,false);x[i][0]=1;}
        var median=QuantileRegressionInference.fitExact(y,x,.5,density);
        assertEquals(Math.PI/(2*n),median.coefficientCovariance()[0],1e-14);
        var d=ExactQuantileTest.data("regular");double[] f=new double[d.y().length];Arrays.fill(f,.3);
        var base=QuantileRegressionInference.fitExact(d.y(),d.x(),.5,f);
        double[] scale={1e9,1e-8,1e4},ys=d.y().clone();
        double[][] xs=Arrays.stream(d.x()).map(double[]::clone).toArray(double[][]::new);
        for(int i=0;i<xs.length;i++){ys[i]*=3;f[i]/=3;for(int j=0;j<3;j++)xs[i][j]*=scale[j];}
        var changed=QuantileRegressionInference.fitExact(ys,xs,.5,f);
        for(int j=0;j<3;j++)for(int k=0;k<3;k++)
            assertEquals(base.coefficientCovariance()[j*3+k],changed.coefficientCovariance()[j*3+k]*scale[j]*scale[k]/9,1e-11);
        double[] copy=base.coefficientCovariance();copy[0]=-1;assertTrue(base.coefficientCovariance()[0]>0);
        copy=base.conditionalDensities();copy[0]=-1;assertEquals(.3,base.conditionalDensities()[0]);
    }
    @Test void invalidDensitiesBandwidthRankAndUncertifiedFitAreRejected() throws Exception {
        var d=ExactQuantileTest.data("regular");double[] f=new double[d.y().length];Arrays.fill(f,.3);
        assertThrows(IllegalStateException.class,()->QuantileRegressionInference.fitExact(d.y(),d.x(),.5,f,new QuantileLinearProgram.Options(1,1e-9)));
        assertThrows(IllegalArgumentException.class,()->QuantileRegressionInference.fitExactIidKernel(d.y(),d.x(),.5,Double.NaN));
        assertThrows(IllegalArgumentException.class,()->QuantileRegressionInference.fitExact(d.y(),d.x(),.5,new double[2]));
        f[0]=0;assertThrows(IllegalArgumentException.class,()->QuantileRegressionInference.fitExact(d.y(),d.x(),.5,f));
        assertThrows(IllegalArgumentException.class,()->QuantileRegressionInference.fitExact(new double[]{1,2,3},new double[][]{{1,1},{1,1},{1,1}},.5,new double[]{.2,.2,.2}));
    }
    private static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
