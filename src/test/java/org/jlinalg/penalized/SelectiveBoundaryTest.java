/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.nio.file.Files;
import java.nio.file.Path;
import jdistlib.Normal;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

final class SelectiveBoundaryTest {
    private static PolyhedralSelectiveInference.Effect constant(double yValue,double xValue,double lambda,double alpha) {
        double[] y=new double[5];Arrays.fill(y,yValue);
        double[][] x=new double[5][1];for(double[] row:x)row[0]=xValue;
        return PolyhedralSelectiveInference.fit(y,x,lambda,alpha,1,false,.95,BackendPolicy.CPU).effects().get(0);
    }

    @Test void nearBoundaryIntervalsMatchIndependentExponentialTailLimit() {
        for(double gap:new double[]{1e-8,1e-10,1e-12}) {
            var e=constant(.4+gap,1,.4,1);
            assertTrue(Double.isFinite(e.confidenceLower()));assertTrue(e.confidenceLower()<e.confidenceUpper());
            // Mills' Q(a)=phi(a)/a*(1+O(a^-2)); here even the smaller a exceeds
            // one million, so this independent limiting inversion is accurate
            // to much better than the 1e-8 relative assertion tolerance.
            double variance=e.untruncatedStandardError()*e.untruncatedStandardError();
            double distance=e.estimate()-e.truncationLower();
            double expectedLow=e.truncationLower()+variance*Math.log(.025)/distance;
            double expectedHigh=e.truncationLower()+variance*Math.log(.975)/distance;
            assertEquals(expectedLow,e.confidenceLower(),Math.abs(expectedLow)*1e-8);
            assertEquals(expectedHigh,e.confidenceUpper(),Math.abs(expectedHigh)*1e-8);
            assertEquals(.975,PolyhedralSelectiveInference.truncatedCdf(e.estimate(),e.confidenceLower(),e.untruncatedStandardError(),e.truncationLower(),e.truncationUpper()),1e-10);
            assertEquals(.025,PolyhedralSelectiveInference.truncatedCdf(e.estimate(),e.confidenceUpper(),e.untruncatedStandardError(),e.truncationLower(),e.truncationUpper()),1e-10);
            var negative=constant(-(.4+gap),1,.4,1);
            assertEquals(-e.confidenceUpper(),negative.confidenceLower(),Math.abs(e.confidenceUpper())*1e-10);
            assertEquals(-e.confidenceLower(),negative.confidenceUpper(),Math.abs(e.confidenceLower())*1e-10);
        }
    }

    @Test void boundaryGapsSurviveHugeMeanShiftsAndFiniteTruncation() {
        assertEquals(-Math.expm1(-1),PolyhedralSelectiveInference.truncatedCdf(1e-12,-1e12,1,0,Double.POSITIVE_INFINITY),1e-14);
        double expected=-Math.expm1(-1)/-Math.expm1(-5);
        assertEquals(expected,PolyhedralSelectiveInference.truncatedCdf(1e-12,-1e12,1,0,5e-12),1e-14);
        assertEquals(1-expected,PolyhedralSelectiveInference.truncatedCdf(-1e-12,1e12,1,-5e-12,0),1e-14);
        assertEquals(.75,PolyhedralSelectiveInference.truncatedCdf(5e-15,0,1,-1e-14,1e-14),1e-14);
    }

    @Test void tailRatioSwitchIsContinuousAndAgreesWithRNormalTails() {
        for(double a:new double[]{0,1,9.999999,10,10.000001,20,100}) {
            for(double d:new double[]{1e-5,.01,1}) {
                double expected=-Math.expm1(Normal.cumulative(-(a+d),0,1,true,true)-Normal.cumulative(-a,0,1,true,true));
                assertEquals(expected,PolyhedralSelectiveInference.truncatedCdf(d,-a,1,0,Double.POSITIVE_INFINITY),2e-12);
            }
        }
        for(double a:new double[]{9.99999,10,10.00001,100}) {
            double hazard=Math.exp(-.5*a*a-.5*Math.log(2*Math.PI)-Normal.cumulative(-a,0,1,true,true));
            double probability=PolyhedralSelectiveInference.truncatedCdf(1e-15,-a,1,0,Double.POSITIVE_INFINITY);
            assertEquals(hazard,probability/1e-15,hazard*1e-11);
        }
    }

    @Test void tinyConstraintProjectionRetainsTheActualSelectionRegion() {
        var baseline=constant(1,1,.4,1);
        for(double scale:new double[]{1e-4,1e-8,1e-16}) {
            // For a single positive column, the event depends only on the
            // L1 threshold. A dominating ridge term cannot remove it.
            var e=constant(1,scale,.4,scale);
            assertEquals(.4,e.truncationLower()*scale,1e-14);
            assertEquals(baseline.pValue(),e.pValue(),1e-13);
            assertEquals(baseline.confidenceLower(),e.confidenceLower()*scale,1e-10);
            assertEquals(baseline.confidenceUpper(),e.confidenceUpper()*scale,1e-10);
            assertEquals(.13660884704067935,e.pValue(),1e-13);
        }
    }

    @Test void positiveAndNegativeExtremeTailsRemainNonzeroAndSymmetric() {
        for(double y:new double[]{4,8,12}) {
            var positive=constant(y,1,.4,1);var negative=constant(-y,1,.4,1);
            double sd=1/Math.sqrt(5);
            double expected=2*Math.exp(Normal.cumulative(-y/sd,0,1,true,true)-Normal.cumulative(-.4/sd,0,1,true,true));
            assertTrue(positive.pValue()>0);assertEquals(expected,positive.pValue(),expected*1e-12);
            assertEquals(positive.pValue(),negative.pValue(),expected*1e-12);
        }
    }

    @Test void unrepresentableVarianceFailsInsteadOfReturningFalseEndpoints() {
        double[] y={1,1,1,1,1};double[][] x={{1},{1},{1},{1},{1}};
        assertThrows(IllegalStateException.class,()->PolyhedralSelectiveInference.fit(y,x,.4,1,
            Double.MIN_VALUE,false,.95,BackendPolicy.CPU));
    }

    @Test void confidenceNearestOneDoesNotRoundItsUpperTailTargetToOne() {
        double[] y={1,1,1,1,1};double[][] x={{1},{1},{1},{1},{1}};
        var e=PolyhedralSelectiveInference.fit(y,x,.4,1,1,false,Math.nextDown(1.0),BackendPolicy.CPU).effects().get(0);
        assertTrue(Double.isFinite(e.confidenceLower()));assertTrue(Double.isFinite(e.confidenceUpper()));
        assertTrue(e.confidenceLower()<e.confidenceUpper());
        double tail=(1-Math.nextDown(1.0))/2;
        assertEquals(tail,PolyhedralSelectiveInference.truncatedCdf(e.estimate(),e.confidenceUpper(),
            e.untruncatedStandardError(),e.truncationLower(),e.truncationUpper()),tail*1e-8);
    }

    @Test void nonorthogonalFixtureMatchesRAndIsInvariantToPredictorUnits() throws Exception {
        var reference=new HashMap<String,List<Double>>();
        for(String line:Files.readAllLines(Path.of("src/test/resources/regression/r-accuracy.csv")).stream().skip(1).toList()) {
            String[] f=line.split(",");reference.computeIfAbsent(f[0],ignored->new ArrayList<>()).add(Double.valueOf(f[2]));
        }
        int n=160;double[] y=new double[n];double[][] x=new double[n][4];
        for(double scale:new double[]{1,1e-8,1e8}) {
            for(int i=0;i<n;i++) {
                double z=-3+6.0*i/(n-1);
                x[i]=new double[]{Math.cos(i*.7)+.2*z,Math.sin(i*.23)-.1*z,Math.cos(i*.33),Math.sin(i*.41)};
                y[i]=.5+.9*x[i][0]-.6*x[i][1]+.5*Math.sin(i*1.13);
                for(int j=0;j<4;j++)x[i][j]*=scale;
            }
            var fit=PolyhedralSelectiveInference.fit(y,x,.12*scale,1,.5,true,.95,BackendPolicy.CPU);
            assertEquals(2,fit.effects().size());
            for(int j=0;j<2;j++) {
                var e=fit.effects().get(j);assertEquals(j,e.predictorIndex());
                assertEquals(reference.get("selective_beta").get(j),e.estimate()*scale,1e-9);
                assertEquals(reference.get("selective_ci").get(2*j),e.confidenceLower()*scale,1e-6);
                assertEquals(reference.get("selective_ci").get(2*j+1),e.confidenceUpper()*scale,1e-6);
                assertEquals(reference.get("selective_limits").get(2*j),e.truncationLower()*scale,1e-6);
                assertEquals(reference.get("selective_limits").get(2*j+1),e.truncationUpper()*scale,1e-6);
            }
        }
    }
}
