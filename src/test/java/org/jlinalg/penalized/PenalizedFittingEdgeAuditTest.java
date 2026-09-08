package org.jlinalg.penalized;

import static org.junit.jupiter.api.Assertions.*;
import java.util.Arrays;
import org.junit.jupiter.api.Test;

class PenalizedFittingEdgeAuditTest {
    @Test void standardizedFitsAreInvariantToTinyPredictorUnitsAndHugeCommonWeights() {
        double[] y = {1,2,5,7,9};
        double[][] x = {{0},{1},{2},{3},{4}};
        PenalizedRegressionResult baseline = PenalizedRegression.lasso(y, x, .1);
        for (double[] weights : new double[][] {{1,1,1,1,1},{1e308,1e308,1e308,1e308,1e308}}) {
            var fit = PenalizedRegression.fit(y, new double[][] {{0},{1e-16},{2e-16},{3e-16},{4e-16}}, .1,
                ElasticNetOptions.builder().alpha(1).observationWeights(weights).build());
            assertTrue(fit.converged());
            assertArrayEquals(baseline.fittedValues(), fit.fittedValues(), 1e-12);
            assertEquals(baseline.objective(), fit.objective(), 1e-12);
        }
    }

    @Test void nearlyExactFitsDoNotLoseResidualSumOfSquares() {
        int n = 40; double[] y = new double[n]; double[][] x = new double[n][2];
        for (int i=0; i<n; i++) {
            x[i][0]=i; x[i][1]=i+.1*Math.sin(i); y[i]=1e8*i+(i%2==0?1:-1);
        }
        var fit = PenalizedRegression.fit(y, x, 0,
            ElasticNetOptions.builder().relativeTolerance(1e-14).build());
        assertTrue(fit.converged());
        double rss = Arrays.stream(fit.residuals()).map(v -> v*v).sum();
        assertTrue(rss > 30);
        assertEquals(rss, fit.weightedResidualSumOfSquares(), 1e-6);
    }

    @Test void cvRejectsIncompletePathsAndPreservesCommonWeightScaling() {
        int n=30; double[] y=new double[n]; double[][] x=new double[n][2];
        for(int i=0;i<n;i++){x[i][0]=i; x[i][1]=i+.3*Math.sin(i); y[i]=2*i+Math.cos(i);}
        assertThrows(IllegalArgumentException.class, () -> PenalizedRegressionCrossValidation.fit(y,x,
            new double[] {1,.1},5,17,ElasticNetOptions.builder().maximumIterations(1).build()));
        double[] w = new double[n]; Arrays.fill(w,1e308);
        var a = PenalizedRegressionCrossValidation.fit(y,x,new double[]{1,.1},5,17,ElasticNetOptions.defaults());
        var b = PenalizedRegressionCrossValidation.fit(y,x,new double[]{1,.1},5,17,
            ElasticNetOptions.builder().observationWeights(w).build());
        assertArrayEquals(a.meanSquaredErrors(),b.meanSquaredErrors(),1e-12);
        assertArrayEquals(a.standardErrors(),b.standardErrors(),1e-12);
    }

    @Test void extremeWeightRatiosKeepRawTrainingWeightsForEachFold() {
        double[] y={0,1,2,3,4,5};double[][] x={{0},{1},{0},{1},{0},{1}};
        var options=ElasticNetOptions.builder().alpha(1)
            .observationWeights(1e200,1e-200,1e-200,1e-200,1e-200,1e-200).build();
        var cv=PenalizedRegressionCrossValidation.fit(y,x,new double[]{100},3,42,options);
        // The dominant held-out row is predicted by the remaining training
        // rows' unweighted mean 3, so the limiting weighted risk is 3^2.
        assertArrayEquals(new double[]{9},cv.meanSquaredErrors(),1e-12);
        assertTrue(Double.isFinite(cv.standardErrors()[0]));
    }
}
