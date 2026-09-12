/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import org.junit.jupiter.api.Test;
final class NonparametricExtensionsTest {
    @Test void automaticLocalInferenceMatchesIndependentWeightedMatrixCalculation() throws Exception {
        var lines=java.nio.file.Files.readAllLines(java.nio.file.Path.of("src/test/resources/estimator-extensions/kernel-data.tsv"));
        double[][] x=new double[lines.size()][2];double[] y=new double[lines.size()];
        for(int i=0;i<y.length;i++){double[] row=Arrays.stream(lines.get(i).split("\\s+")).mapToDouble(Double::parseDouble).toArray();x[i][0]=row[0];x[i][1]=row[1];y[i]=row[2];}
        double[] expected=java.nio.file.Files.readAllLines(java.nio.file.Path.of("src/test/resources/estimator-extensions/kernel-reference.tsv")).stream().mapToDouble(Double::parseDouble).toArray();
        var types=new ProductKernelRegression.Type[]{ProductKernelRegression.Type.CONTINUOUS,ProductKernelRegression.Type.UNORDERED};
        var fit=ProductKernelRegression.infer(x,y,new double[]{.5,1},types,.95);
        assertEquals(expected[0],fit.estimate(),1e-11);assertEquals(expected[1],fit.standardError(),1e-11);assertEquals(expected[2],fit.bandwidths()[0],1e-13);
    }
    @Test void mixedLocalLinearReproducesPlaneWithinCategory() {
        double[][] x=new double[60][3];double[] y=new double[60];
        for(int i=0;i<60;i++){x[i]=new double[]{i%5,i/5%6,i%2};y[i]=2+3*x[i][0]-2*x[i][1]+8*x[i][2];}
        var types=new ProductKernelRegression.Type[]{ProductKernelRegression.Type.CONTINUOUS,ProductKernelRegression.Type.CONTINUOUS,ProductKernelRegression.Type.UNORDERED};
        assertEquals(2+3*1.2-2*2.7+8,ProductKernelRegression.predict(x,y,new double[]{1.2,2.7,1},types,new double[]{.8,1.2,0}),1e-11);
        assertThrows(IllegalArgumentException.class,()->ProductKernelRegression.predict(x,y,new double[]{1,2,3},types,new double[]{1,1,0}));
    }
    @Test void automaticBandwidthAndInferenceRespectResponseAffineTransformation() {
        int n=160;double[][] x=new double[n][2];double[] y=new double[n],scaled=new double[n];
        for(int i=0;i<n;i++){x[i]=new double[]{(i%80)/79.0,i/80};y[i]=2+x[i][0]+Math.sin(13*i)*.3;scaled[i]=10+4*y[i];}
        var types=new ProductKernelRegression.Type[]{ProductKernelRegression.Type.CONTINUOUS,ProductKernelRegression.Type.ORDERED};
        var a=ProductKernelRegression.infer(x,y,new double[]{.5,0},types,.95);
        var b=ProductKernelRegression.infer(x,scaled,new double[]{.5,0},types,.95);
        assertEquals(10+4*a.estimate(),b.estimate(),1e-10);assertEquals(4*a.standardError(),b.standardError(),1e-10);
        assertArrayEquals(a.bandwidths(),b.bandwidths());assertEquals(80,a.stratumSize());assertTrue(a.standardError()>0);
    }
    @Test void massPointIntervalsHaveConservativeExactCoverage() {
        // Enumerate Binomial(20,.7) counts: population median is one.
        double coverage=0;
        for(int ones=0;ones<=20;ones++) {
            double[] y=new double[20];Arrays.fill(y,0,ones,1);
            var interval=QuantileRegressionInference.marginalInterval(y,.5,.95);
            if(interval.lower()<=1&&interval.upper()>=1)coverage+=jdistlib.Binomial.density(ones,20,.7,false);
        }
        assertTrue(coverage>=.95);assertTrue(coverage<=1+1e-12);
    }
    @Test void automaticQuantileDensityAndClusterCovarianceAreFinite() {
        int n=90;double[] y=new double[n],density=new double[n];double[][] x=new double[n][1];int[] groups=new int[n];
        for(int i=0;i<n;i++){y[i]=Math.sin(i*1.713)+.001*i;x[i][0]=1;groups[i]=i/3;density[i]=.4;}
        var auto=QuantileRegressionInference.fitExactIidKernel(y,x,.5);
        var cluster=QuantileRegressionInference.fitExactCluster(y,x,.5,density,groups);
        var hac=QuantileRegressionInference.fitExactHac(y,x,.5,density,4);
        double[] psi=new double[n];double mean=0,meat=0;
        for(int i=0;i<n;i++){psi[i]=.5-(hac.fit().residuals()[i]<0?1:0);mean+=psi[i]/n;}
        for(int i=0;i<n;i++){psi[i]-=mean;meat+=psi[i]*psi[i];}
        for(int lag=1;lag<=4;lag++)for(int i=lag;i<n;i++)meat+=2*(1-lag/5.0)*psi[i]*psi[i-lag];
        assertEquals(meat/Math.pow(n*.4,2),hac.coefficientCovariance()[0],1e-12);
        assertTrue(auto.bandwidth()>0);assertTrue(auto.standardErrors()[0]>0);assertTrue(cluster.standardErrors()[0]>0);
    }
}
