/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class PartiallyLinearNoiseTest {
    @Test void periodicSmootherGuardsItsMinimumWindowSize() {
        assertThrows(IllegalArgumentException.class,()->SuperSmoother.fit(
            new double[]{0,.3,.6,1},new double[]{1,2,3,4},null,0,true,0));
        double[] x=new double[5],y=new double[5];
        for(int i=0;i<5;i++){x[i]=i/4.0;y[i]=Math.cos((i+1)*1.7);}
        assertArrayEquals(new double[]{-.1299327931310073,-.1385516120479152,
            -.1318166376041629,-.1321979937678831,-.1465499486326298},
            SuperSmoother.fit(x,y,null,0,true,0).fittedValues(),1e-12);
    }
    @Test void kernelWeightedMeanDoesNotOverflowAnUnnormalizedSum() {
        double value=Double.MAX_VALUE/2;
        assertEquals(value,KernelRegression.predict(new double[]{0,0,0,0},
            new double[]{value,value,value,value},new double[]{0},1)[0]);
    }
    @Test void smoothingInducedNoiseCovarianceIsNotTreatedAsIndependent() {
        // Each pair is smoothed to its mean. The slope uses the independent
        // original pair contrast, whose variance is exactly 1/n at unit noise.
        for(int n:new int[]{8,80}) {
            double[] z=new double[n];double[][] x=new double[n][1];
            for(int i=0;i<n;i++){z[i]=100*(i/2);x[i][0]=i%2==0?1:-1;}
            double expectedHc3=0,expectedHomoskedastic=0,exactVariance=0;
            for(int j=0;j<n;j++) {
                double[] y=new double[n];y[j]=1;
                var fit=PartiallyLinearInference.fit(y,x,z,1,BackendPolicy.CPU);
                expectedHc3+=fit.slopeCovariance()[0];
                expectedHomoskedastic+=fit.homoskedasticSlopeCovariance()[0];
                exactVariance+=Math.pow(fit.fit().coefficients()[0],2);
                assertEquals(n/2.0-1,fit.residualNoiseDegreesOfFreedom(),1e-10);
            }
            assertEquals(1.0/n,exactVariance,1e-12);
            assertEquals(exactVariance,expectedHomoskedastic,1e-12);
            // Delete-one generalized HC3 is conservative for these high-leverage
            // pair smoothers, not an exact finite-sample variance estimator.
            assertEquals(2.0/(n-2),expectedHc3,1e-12);
            assertTrue(expectedHc3>=exactVariance);
        }
    }
}
