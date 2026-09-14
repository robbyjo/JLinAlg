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
        assertTrue(fit.jointParameterInferenceAvailable());
        assertArrayEquals(new double[]{fit.coefficients()[0],fit.coefficients()[1],
            fit.effectiveAutoregressive()[0]},fit.jointParameterEstimates(),1e-12);
        double[] expectedJoint={
            .014446134660396506,-.001908199556996274,.000901338395385856,
            -.001908199556996274,.026580742091447916,-.000820899733570494,
            .000901338395385856,-.000820899733570494,.012504040728895988};
        assertArrayEquals(expectedJoint,fit.jointParameterCovariance(),4e-4);
        var conditional=fit.forecastConditional(new double[][]{{1,.2},{1,.4}},.95);
        var uncertain=fit.forecastWithParameterUncertainty(new double[][]{{1,.2},{1,.4}},.95);
        assertArrayEquals(conditional.means(),uncertain.means(),1e-12);
        for(int i=0;i<2;i++)assertTrue(uncertain.standardErrors()[i]>conditional.standardErrors()[i]);
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
        assertTrue(fit.jointParameterInferenceAvailable());
        assertArrayEquals(fit.conditionalCoefficientCovariance(),fit.jointParameterCovariance(),2e-6);
        var conditional=fit.forecastConditional(new double[][]{{1,6}},.95);
        var uncertain=fit.forecastWithParameterUncertainty(new double[][]{{1,6}},.95);
        double[] c=fit.jointParameterCovariance();double parameterVariance=c[0]+12*c[1]+36*c[3];
        assertEquals(conditional.standardErrors()[0]*conditional.standardErrors()[0]+parameterVariance,
            uncertain.standardErrors()[0]*uncertain.standardErrors()[0],2e-6);
        assertThrows(IllegalArgumentException.class,()->fit.forecast(new double[][]{{6}},.95));
    }
    @Test void conditionalOperationsDoNotComputeJointHessian() {
        double[] y={2,1,4,3,6,5};double[][] x={{1,0},{1,1},{1,2},{1,3},{1,4},{1,5}};
        var fit=ArimaRegression.fit(y,x,new ArimaOrder(0,0,0),null);
        assertFalse(fit.jointParameterInferenceComputed());
        fit.forecastConditional(new double[][]{{1,6}},.95);
        fit.smoothErrors();
        assertEquals(2,fit.jointParameterEstimates().length);
        assertFalse(fit.jointParameterInferenceComputed());
        assertTrue(fit.jointParameterInferenceAvailable());
        assertTrue(fit.jointParameterInferenceComputed());
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
    @Test void covarianceToleranceDoesNotSquareLargeInnovationScale() {
        double scale=1e12;
        double[] roundoff={-5e-9*scale};
        ArimaSmoothing.stabilizeCovariance(roundoff,1,scale);
        assertEquals(0,roundoff[0],0);
        assertThrows(IllegalArgumentException.class,()->
            ArimaSmoothing.stabilizeCovariance(
                new double[]{-1e-4*scale},1,scale));
    }
    @Test void differencingRemovesUnidentifiedIntercept() {
        double[] y={1,3,2,6,4,7};double[][] x={{1},{1},{1},{1},{1},{1}};
        assertThrows(IllegalArgumentException.class,()->ArimaRegression.fit(y,x,new ArimaOrder(0,1,0),null));
    }
    @Test void exactDiffuseBackwardSmootherScalesBeyondFormerDateBound() {
        int dates=5_000;double[] y=new double[dates];java.util.Arrays.fill(y,Double.NaN);
        y[0]=2;y[dates-1]=8;
        var result=ArimaSmoothing.smooth(y,new double[0],new double[0],new double[]{1,-1},2);
        int t=1_937;double fraction=(double)t/(dates-1);
        assertEquals(2+6*fraction,result.signal()[t],2e-10);
        assertEquals(2.0*t*(dates-1-t)/(dates-1),result.signalVariance()[t],2e-7);
        assertEquals(0,result.signalVariance()[0],1e-12);
        assertEquals(0,result.signalVariance()[dates-1],1e-12);
        assertEquals(dates,result.states().length);
    }
    @Test void jointEstimatesExposeSeasonalCoordinatesInsteadOfEffectivePolynomial() {
        java.util.Random random=new java.util.Random(20260914L);
        double[] full=new double[520];
        for(int t=0;t<full.length;t++)
            full[t]=random.nextGaussian()+(t>=4?.55*full[t-4]:0);
        int n=320;double[] y=new double[n];double[][] x=new double[n][2];
        for(int i=0;i<n;i++) {
            int t=i+200;y[i]=2+.01*t+full[t];x[i]=new double[]{1,t};
        }
        var options=ArimaOptions.builder()
            .seasonalOrder(SeasonalArimaOrder.of(1,0,0,4))
            .includeMean(false).build();
        var fit=ArimaRegression.fit(y,x,ArimaOrder.arma(0,0),options);
        assertTrue(fit.converged());
        assertTrue(fit.jointParameterInferenceAvailable());
        double[] estimates=fit.jointParameterEstimates();
        assertEquals(3,estimates.length);
        assertArrayEquals(fit.coefficients(),
            java.util.Arrays.copyOf(estimates,2),1e-12);
        assertEquals(.55,estimates[2],.12);
        assertEquals(9,fit.jointParameterCovariance().length);
        assertEquals(4,fit.effectiveAutoregressive().length);
        assertEquals(estimates[2],fit.effectiveAutoregressive()[3],1e-12);
    }
}
