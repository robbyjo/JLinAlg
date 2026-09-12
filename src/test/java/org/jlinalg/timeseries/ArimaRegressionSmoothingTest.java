/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

final class ArimaRegressionSmoothingTest {
    private static double[][] matrix(String name) throws Exception {
        return java.nio.file.Files.readAllLines(java.nio.file.Path.of("src/test/resources/estimator-extensions/"+name)).stream()
            .map(s->java.util.Arrays.stream(s.split("\\s+")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
    }
    @Test void regressionAndHistoricalStatesMatchR() throws Exception {
        double[][] data=matrix("arima-data.tsv");double[] y=new double[data.length];double[][] x=new double[data.length][2];
        for(int i=0;i<y.length;i++){y[i]=data[i][0];x[i][0]=data[i][1];x[i][1]=data[i][2];}
        var fit=ArimaRegression.fit(y,x,new ArimaOrder(1,0,0),null);double[][] reference=matrix("arima-reference.tsv");
        assertTrue(fit.converged());assertEquals(reference[1][0],fit.coefficients()[0],2e-5);assertEquals(reference[2][0],fit.coefficients()[1],2e-5);
        assertEquals(reference[3][0],fit.innovationVariance(),2e-5);assertEquals(reference[4][0],fit.logLikelihood(),1e-7);
        for(boolean integrated:new boolean[]{false,true}) {
            var sm=ArimaSmoothing.smooth(y,new double[]{.4},new double[]{.2},integrated?new double[]{1,-1}:new double[]{1},1);
            double[][] expected=matrix(integrated?"smooth-integrated.tsv":"smooth-stationary.tsv");
            for(int i=0;i<y.length;i++)assertArrayEquals(expected[i],sm.states()[i],integrated?2e-5:1e-9);
        }
    }
    @Test void whiteNoiseRegressionIsOlsAndForecastUsesFutureDesign() {
        double[] y={2,1,4,3,6,5};double[][] x={{1,0},{1,1},{1,2},{1,3},{1,4},{1,5}};
        var fit=ArimaRegression.fit(y,x,new ArimaOrder(0,0,0),null);
        assertTrue(fit.converged());assertArrayEquals(new double[]{10.0/7,29.0/35},fit.coefficients(),1e-10);
        assertEquals(10.0/7+6*29.0/35,fit.forecast(new double[][]{{1,6}},.95).means()[0],1e-10);
        assertThrows(IllegalArgumentException.class,()->fit.forecast(new double[][]{{6}},.95));
    }
    @Test void randomWalkBridgeUsesBothEndpointsAndDiffusePrefix() {
        var result=ArimaSmoothing.smooth(new double[]{Double.NaN,2,Double.NaN,Double.NaN,8,Double.NaN},
            new double[0],new double[0],new double[]{1,-1},2);
        assertArrayEquals(new double[]{2,2,4,6,8,8},result.signal(),1e-9);
        assertArrayEquals(new double[]{2,0,4.0/3,4.0/3,0,2},result.signalVariance(),1e-9);
        result.states()[0][0]=999;assertNotEquals(999,result.states()[0][0]);
    }
    @Test void stationaryArBridgeMatchesAnalyticConditionalNormal() {
        double phi=.5;
        var fit=ArimaSmoothing.smooth(new double[]{1,Double.NaN,3},new double[]{phi},new double[0],new double[]{1},2);
        assertEquals(phi*(1+3)/(1+phi*phi),fit.signal()[1],1e-12);
        assertEquals(2/(1+phi*phi),fit.signalVariance()[1],1e-12);
    }
    @Test void differencingRemovesUnidentifiedIntercept() {
        double[] y={1,3,2,6,4,7};double[][] x={{1},{1},{1},{1},{1},{1}};
        assertThrows(IllegalArgumentException.class,()->ArimaRegression.fit(y,x,new ArimaOrder(0,1,0),null));
    }
}
