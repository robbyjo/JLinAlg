package org.jlinalg.ols;

import static org.junit.jupiter.api.Assertions.*;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

class OlsEdgeAuditTest {
    @Test void predictorUnitsDoNotDetermineFullRankOrInference() {
        double[] y = {1, 2, 5, 7, 9};
        for (double scale : new double[] {1e-16, 1, 1e16}) {
            double[][] x = {{1, 0}, {1, scale}, {1, 2*scale}, {1, 3*scale}, {1, 4*scale}};
            OlsResult fit = Ols.fit(y, x, OlsOptions.defaults(), BackendPolicy.CPU);
            assertEquals(2, fit.rank());
            assertEquals(2.1, fit.coefficients()[1] * scale, 1e-12);
            assertEquals(Math.sqrt(.7/30), fit.standardErrors()[1] * scale, 1e-12);
            assertEquals(0.0008328504718253874, fit.pValues()[1], 1e-14);
        }
    }

    @Test void originalCoordinateMinimumNormAndEstimabilityArePreserved() {
        double[][] x = {{1, 0, 0}, {1, 1, 2}, {1, 2, 4}, {1, 3, 6}, {1, 4, 8}};
        OlsResult fit = Ols.fit(new double[] {1, 2, 5, 7, 9}, x,
            new OlsOptions(RankDeficiencyStrategy.MINIMUM_NORM, .95), BackendPolicy.CPU);
        assertArrayEquals(new double[] {.6, .42, .84}, fit.coefficients(), 1e-12);
        assertTrue(Double.isFinite(fit.pValues()[0]));
        assertTrue(Double.isNaN(fit.standardErrors()[1]));
        assertTrue(Double.isNaN(fit.associationStatistics().pValues()[1]));
        assertThrows(IllegalArgumentException.class, () -> fit.testContrast(new double[][] {{0, 1, 0}}));
        assertEquals(0.0008328504718253874,
            fit.testContrast(new double[][] {{0, 1, 2}}).pValue(), 1e-12);
        // The mean-norm covariance remains available for estimable contrasts.
        double[] c = fit.covariance();
        assertEquals(.7/30, c[4] + 4*c[5] + 4*c[8], 1e-12);
    }

    @Test void fastPathBoundaryChecksColumnRatioNotOnlyIndividualExtremes() {
        int n=10000; double[][] x=new double[n][2], reference=new double[n][2];double[] y=new double[n];
        for(int i=0;i<n;i++) {
            double z=i%2==0?1:-1;x[i][0]=1e-6;x[i][1]=1e6*z;
            reference[i][0]=1;reference[i][1]=z;y[i]=2+3*z+.1*Math.sin(i);
        }
        var expected=Ols.fit(y,reference,OlsOptions.defaults(),BackendPolicy.CPU);
        var fit=Ols.fit(y,x,OlsOptions.defaults(),BackendPolicy.CPU);
        assertEquals(2,fit.rank());
        assertEquals(expected.coefficients()[0],fit.coefficients()[0]*1e-6,1e-12);
        assertEquals(expected.coefficients()[1],fit.coefficients()[1]*1e6,1e-12);
        assertArrayEquals(expected.fittedValues(),fit.fittedValues(),1e-12);
    }

    @Test void wideOriginalCoordinateMinimumNormIsPreserved() {
        double[][] x={{1,2,0,0,0},{1,2,0,0,0},{1,2,0,0,0},{1,2,0,0,0}};
        var fit=Ols.fit(new double[]{0,0,4,8},x,
            new OlsOptions(RankDeficiencyStrategy.MINIMUM_NORM,.95),BackendPolicy.CPU);
        assertEquals(1,fit.rank());assertEquals(44,fit.residualSumOfSquares(),1e-10);
        assertArrayEquals(new double[]{.6,1.2,0,0,0},fit.coefficients(),1e-10);
    }

    @Test void inaccurateOriginalSvdCannotSilentlyDropOrBiasAnEstimableIntercept() {
        for(int columns:new int[]{3,5}) {
            double[][] x=new double[4][columns];double[] z={-2,-1,1,2};
            for(int i=0;i<4;i++){x[i][0]=1;for(int j=1;j<columns;j++)x[i][j]=z[i]*1e16;}
            var error=assertThrows(IllegalArgumentException.class,()->Ols.fit(new double[]{0,0,4,8},x,
                new OlsOptions(RankDeficiencyStrategy.MINIMUM_NORM,.95),BackendPolicy.CPU));
            assertTrue(error.getMessage().contains("original-coordinate"),error.getMessage());
        }
    }
}
