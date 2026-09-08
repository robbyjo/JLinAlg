package org.jlinalg.gee;

import static org.junit.jupiter.api.Assertions.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.junit.jupiter.api.Test;

class GeeScaleAuditTest {
    @Test void gammaUnitsDoNotChangeEstimatingEquationsOrSandwich() {
        for (GeeCorrelation correlation : new GeeCorrelation[] {GeeCorrelation.INDEPENDENCE, GeeCorrelation.EXCHANGEABLE}) {
            GeeResult a = fit(1, correlation, 1), b = fit(1e-8, correlation, 1);
            assertTrue(a.converged(), a.convergenceMessage());
            assertTrue(b.converged(), b.convergenceMessage());
            assertEquals(a.coefficients()[0] + Math.log(1e-8), b.coefficients()[0], 1e-8);
            assertEquals(a.coefficients()[1], b.coefficients()[1], 1e-8);
            assertEquals(a.dispersion(), b.dispersion(), 1e-8);
            assertArrayEquals(a.covariance(), b.covariance(), 1e-8);
        }
    }

    @Test void serialAndParallelAccumulationsAgree() {
        GeeResult a=fit(1, GeeCorrelation.INDEPENDENCE, 1);
        GeeResult b=fit(1, GeeCorrelation.INDEPENDENCE, 2);
        assertTrue(a.converged()); assertTrue(b.converged());
        assertArrayEquals(a.coefficients(),b.coefficients(),1e-12);
        assertArrayEquals(a.covariance(),b.covariance(),1e-12);
    }

    @Test void gaussianDispersionAndNaiveCovarianceFollowResponseUnits() {
        var fit=Gee.fit(new double[]{1e-12,2e-12,3e-12,4e-12,5e-12,6e-12},
            new double[][]{{1},{1},{1},{1},{1},{1}},new int[]{1,1,2,2,3,3},null,
            GlmFamilies.gaussian(),null,null,GeeOptions.builder().covariance(GeeCovariance.NAIVE).build(),BackendPolicy.CPU);
        assertTrue(fit.converged(),fit.convergenceMessage());
        assertEquals(3.5e-24,fit.dispersion(),1e-36);
        assertEquals(3.5e-24/6,fit.covariance()[0],1e-36);
    }

    private static GeeResult fit(double scale, GeeCorrelation correlation, int parallelism) {
        int n=60; double[] y=new double[n]; double[][] x=new double[n][2]; int[] id=new int[n];
        for(int i=0;i<n;i++){x[i][0]=1;x[i][1]=i%4-1.5;y[i]=scale*Math.exp(.3+.5*x[i][1]+.7*Math.sin(i));id[i]=i/3;}
        return Gee.fit(y,x,id,null,GlmFamilies.gamma(),null,null,
            GeeOptions.builder().correlation(correlation).parallelism(parallelism).parallelThreshold(1).build(),BackendPolicy.CPU);
    }
}
