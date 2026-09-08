/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import org.jlinalg.compute.BackendPolicy;
import org.junit.jupiter.api.Test;

final class SelectionAccuracyTest {
    @Test void fullSampleWeightsAreSplitAndPropagated() {
        double[] y=new double[20],w=new double[20];double[][] x=new double[20][1];
        for(int i=0;i<20;i++){y[i]=i+Math.sin(i);x[i][0]=i;w[i]=1;}
        var fit=SelectionAwarePenalizedInference.fit(y,x,.1,ElasticNetOptions.builder().observationWeights(w).build(),.5,1e-8,BackendPolicy.CPU);
        assertEquals(10,fit.selectionObservations());assertEquals(10,fit.inferenceObservations());
    }
    @Test void highDimensionalSelectionDoesNotRequireFullDesignRank() {
        double[] y=new double[20];double[][] x=new double[20][30];
        for(int i=0;i<20;i++){y[i]=Math.sin(i);for(int j=0;j<30;j++)x[i][j]=Math.cos((i+1)*(j+1));}
        var fit=SelectionAwarePenalizedInference.fit(y,x,100,ElasticNetOptions.defaults(),.5,1e-8,BackendPolicy.CPU);
        assertEquals(10,fit.selectionObservations());assertEquals(0,fit.activePredictorIndices().length);
    }
    @Test void fixedPenaltyConditionalIntervalHasAnalyticTruncation() {
        double[] y={1.8,2.1,1.9,2.3,1.9};double[][] x={{1},{1},{1},{1},{1}};
        for(double alpha:new double[]{1,.5}) {
            var fit=PolyhedralSelectiveInference.fit(y,x,.4,alpha,1,false,.95,BackendPolicy.CPU);
            var effect=fit.effects().get(0);
            assertEquals(2,effect.estimate(),1e-12);
            assertEquals(.4*alpha,effect.truncationLower(),1e-10);
            assertEquals(Double.POSITIVE_INFINITY,effect.truncationUpper());
            assertTrue(effect.confidenceLower()<2);assertTrue(effect.confidenceUpper()>2);
            assertEquals(.975,PolyhedralSelectiveInference.truncatedCdf(2,effect.confidenceLower(),1/Math.sqrt(5),.4*alpha,Double.POSITIVE_INFINITY),1e-8);
            assertEquals(.025,PolyhedralSelectiveInference.truncatedCdf(2,effect.confidenceUpper(),1/Math.sqrt(5),.4*alpha,Double.POSITIVE_INFINITY),1e-8);
        }
    }
    @Test void conditionalPValuesAreCalibratedOnTheTruncatedNormal() {
        double sd=1/Math.sqrt(5),threshold=.4;
        double base=jdistlib.Normal.cumulative(threshold,0,sd,true,false);
        for(double u:new double[]{.001,.025,.2,.5,.8,.975,.999}) {
            double t=jdistlib.Normal.quantile(base+u*(1-base),0,sd,true,false);
            double[] y={t,t,t,t,t};double[][] x={{1},{1},{1},{1},{1}};
            var e=PolyhedralSelectiveInference.fit(y,x,threshold,1,1,false,.95,BackendPolicy.CPU).effects().get(0);
            assertEquals(2*Math.min(u,1-u),e.pValue(),1e-11);
        }
        double[] y={5,5,5,5,5};double[][] x={{1},{1},{1},{1},{1}};
        assertTrue(PolyhedralSelectiveInference.fit(y,x,threshold,1,1,false,.95,BackendPolicy.CPU).effects().get(0).pValue()>0);
    }
}
