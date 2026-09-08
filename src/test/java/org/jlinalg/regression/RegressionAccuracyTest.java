/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.regression;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

final class RegressionAccuracyTest {
    @Test void medianHasCorrectStationaryPoint() {
        var fit = QuantileRegression.fit(new double[]{1,2,3}, new double[][]{{1},{1},{1}}, .5);
        assertEquals(2, fit.coefficients()[0], 1e-6);
        assertTrue(fit.converged());
    }
    @Test void quantileIsInvariantToPredictorUnits() {
        var fit = QuantileRegression.fit(new double[]{1,2,3}, new double[][]{{1e12},{1e12},{1e12}}, .5);
        assertEquals(2, fit.coefficients()[0]*1e12, 1e-6);
        assertTrue(fit.converged());
    }
    @Test void multinomialIsInvariantToPredictorUnits() {
        var fit = MultinomialRegression.fit(new int[]{0,1,1,1}, new double[][]{{1e12},{1e12},{1e12},{1e12}}, 2);
        assertEquals(.75, fit.probabilities()[1], 1e-7);
        assertEquals(Math.log(3), fit.coefficients()[0]*1e12, 1e-6);
        assertTrue(fit.converged());
    }
    @Test void kernelDoesNotUnderflowAtRemoteQuery() {
        assertEquals(2, KernelRegression.fit(new double[]{0,1}, new double[]{1,2}, 1).predict(new double[]{1000})[0], 1e-12);
    }
    @Test void smootherPredictionsAreBatchInvariant() {
        double[] x = new double[40], y = new double[40];
        for(int i=0;i<x.length;i++){x[i]=i/10.;y[i]=Math.sin(x[i])+.3*Math.cos(i*1.7);}
        var fit=SuperSmoother.fit(x,y);
        double expected=fit.predict(new double[]{1.23})[0];
        for(double value:fit.predict(new double[]{1.23,1.23,1.23,1.23})) assertEquals(expected,value,0);
    }
    @Test void smootherTiesDoNotCrash() {
        var fit=SuperSmoother.fit(new double[]{0,0,0,0},new double[]{1,2,3,4});
        assertEquals(2.5, fit.predict(new double[]{10})[0], 1e-12);
    }
}
